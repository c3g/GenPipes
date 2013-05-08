package ca.mcgill.genome;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.io.input.CountingInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import au.com.bytecode.opencsv.CSVReader;
import ca.mcgill.genome.SampleAlleleCounts.ReadBases;

public class VariantAlleleFrequency implements Runnable {
	private static final Logger log = LoggerFactory.getLogger(VariantAlleleFrequency.class);
	private final File output;
	private final File inputCSV;
	private final Map<String, SampleAlleleCounts> samplesBAM;
	private final String chromosome;
	private final byte threads;
	private PairMapping pairs = null;

	public VariantAlleleFrequency(File output, File inputCSV, Map<String, File> samplesFile, String chromosome, byte threads) {
		this.output = output;
		this.inputCSV = inputCSV;
		this.samplesBAM = new HashMap<String, SampleAlleleCounts>();
		this.chromosome = chromosome;
		this.threads = threads;
		
		for(String sampleName : samplesFile.keySet()) {
			samplesBAM.put(sampleName, new SampleAlleleCounts(sampleName, samplesFile.get(sampleName), 15, 0));
		}
	}

	public void setPairs(PairMapping pairs) {
		this.pairs = pairs;
	}

	@Override
	public void run() {
		ExecutorService bamCoverageExecutor = Executors.newFixedThreadPool(threads);

		CSVReader reader = null;
		PrintWriter writer = null;
		try {
			long inputCSVSize = inputCSV.length();
			CountingInputStream countInput = new CountingInputStream(new BufferedInputStream(new FileInputStream(inputCSV), 16*1024*1024));
			reader = new CSVReader(new InputStreamReader(countInput, Charset.forName("ASCII")), '\t');
			writer = new PrintWriter(output, "ASCII");
			String[] nextLine = reader.readNext();
			if(nextLine == null)
				throw new RuntimeException("CSV is empty");

			int chromosomeIdx = -1;
			int positionIdx = -1;
			Map<String, Integer> sampleToColumn = new HashMap<String, Integer>();
			Map<Integer, String> columnToSample = new HashMap<Integer, String>();
			for(int idx=0; idx < nextLine.length; idx++) {
				if(nextLine[idx].equals("chromosome")) {
					chromosomeIdx = idx;
				}
				else if(nextLine[idx].equals("position")) {
					positionIdx = idx;
				}
				else {
					String sample = nextLine[idx];
					if(pairs != null) {
						String normal = pairs.getNormalFromSample(sample);
						if(normal != null) {
							if(samplesBAM.containsKey(normal)) {
								sampleToColumn.put(sample, idx);
								columnToSample.put(idx, sample);
							}
						}
						String tumor = pairs.getTumorFromSample(sample);
						if(tumor != null) {
							if(samplesBAM.containsKey(tumor)) {
								sampleToColumn.put(sample, idx);
								columnToSample.put(idx, sample);
							}
						}
					}
					else {
						if(samplesBAM.containsKey(sample)) {
							sampleToColumn.put(sample, idx);
							columnToSample.put(idx, sample);
						}
					}
				}
			}

			for(int idx=0; idx < nextLine.length; idx++) {
				if(idx > 0) {
					writer.print('\t');
				}
				String sample = columnToSample.get(idx);
				if(sample != null) {
					writer.print(sample + "allele counts");
					writer.print('\t');
				}
				writer.print(nextLine[idx]);
			}
			writer.println();
			
			Map<String, Future<EnumMap<ReadBases, Integer>>> sampleFutures = new HashMap<String, Future<EnumMap<ReadBases, Integer>>>(sampleToColumn.size());
			long prevPercent=0;
			while ((nextLine = reader.readNext()) != null) {
				long percent = countInput.getByteCount()*100l/inputCSVSize;
				if(percent != prevPercent) {
					log.info("Completed: {}%", percent);
					prevPercent = percent;
				}
				sampleFutures.clear();
				String chromosome = nextLine[chromosomeIdx];
				// Skip chromosomes that don't match if asked for.
				if(this.chromosome != null && !this.chromosome.equals(chromosome))
					continue;

				int position = Integer.parseInt(nextLine[positionIdx]);

				for(String sampleName : sampleToColumn.keySet()) {
					if(pairs != null) {
						String normal = pairs.getNormalFromSample(sampleName);
						Future<EnumMap<ReadBases, Integer>> future;
						future = bamCoverageExecutor.submit(samplesBAM.get(normal).getCallableAlleleCountsComputer(chromosome, position));
						sampleFutures.put(normal, future);
						String tumor = pairs.getTumorFromSample(sampleName);
						future = bamCoverageExecutor.submit(samplesBAM.get(tumor).getCallableAlleleCountsComputer(chromosome, position));
						sampleFutures.put(tumor, future);
					}
					else {
						Future<EnumMap<ReadBases, Integer>> future = bamCoverageExecutor.submit(samplesBAM.get(sampleName).getCallableAlleleCountsComputer(chromosome, position));
						sampleFutures.put(sampleName, future);
					}
				}
				
				Map<String, StringBuilder> toPrint = new HashMap<String, StringBuilder>();
				for(String sampleName : sampleFutures.keySet()) {
					Future<EnumMap<ReadBases, Integer>> future = sampleFutures.get(sampleName);
					EnumMap<ReadBases, Integer> baseCounts = future.get();
					
					if(pairs != null) {
						if(pairs.getSampleFromNormal(sampleName) != null) {
							String sample = pairs.getSampleFromNormal(sampleName);
							if(!toPrint.containsKey(sample))
								toPrint.put(sample, new StringBuilder());

							StringBuilder sb = new StringBuilder();
							boolean first = true;
							for(ReadBases base : baseCounts.keySet()) {
								if(first) {
									first = false;
								}
								else {
									sb.append(' ');
								}
								sb.append(base.toString());
								sb.append(':');
								sb.append(baseCounts.get(base));
							}
							sb.append(" , ");
							toPrint.get(sample).insert(0, sb.toString());
						}
						else if(pairs.getSampleFromTumor(sampleName) != null) {
							String sample = pairs.getSampleFromTumor(sampleName);
							if(!toPrint.containsKey(sample))
								toPrint.put(sample, new StringBuilder());

							StringBuilder sb = toPrint.get(sample);
							boolean first = true;
							for(ReadBases base : baseCounts.keySet()) {
								if(first) {
									first = false;
								}
								else {
									sb.append(' ');
								}
								sb.append(base.toString());
								sb.append(':');
								sb.append(baseCounts.get(base));
							}
						}
						else {
							throw new RuntimeException("Unknown sample: " + sampleName);
						}
					}
					else {
						StringBuilder sb = new StringBuilder();
						boolean first = true;
						for(ReadBases base : baseCounts.keySet()) {
							if(first) {
								first = false;
							}
							else {
								sb.append(' ');
							}
							sb.append(base.toString());
							sb.append(':');
							sb.append(baseCounts.get(base));
						}
						toPrint.put(sampleName, sb);
					}
					
				} // for futures
				
				for(int idx=0; idx < nextLine.length; idx++) {
					if(idx > 0) {
						writer.print('\t');
					}
					String sample = columnToSample.get(idx);
					if(sample != null) {
						writer.print(toPrint.get(sample));
						writer.print('\t');
					}
					writer.print(nextLine[idx]);
				}
				writer.println();
			} // while line reader
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} catch (ExecutionException e) {
			throw new RuntimeException(e);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
			if(writer != null) {
				writer.close();
			}
			
			bamCoverageExecutor.shutdownNow();
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Options options = new Options();
		options.addOption("h", "help", false, "print this message");
		options.addOption("o", "output", true, "output csv");
		options.addOption("i", "input", true, "input csv");
		options.addOption("p", "pairs", true, "pair file");
		options.addOption("s", "samplesFile", true, "input file containing sample names and bams");
		options.addOption("t", "threads", true, "nbThreads");
		options.addOption("r", "region", true, "region (chromosome) to annoatate");

		CommandLineParser parser = new PosixParser();
		try {
			CommandLine line = parser.parse(options, args);
			long start = System.currentTimeMillis();
			if (line.hasOption("help") || !line.hasOption("output") || !line.hasOption("input")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp(VariantAlleleFrequency.class.getName(), options);
			} else {
				byte threads = Byte.parseByte(line.getOptionValue("threads"));
				File output = new File(line.getOptionValue("output"));
				File inputCSV = new File(line.getOptionValue("input"));
				Map<String, File> samples = parseSamplesFile(new File(line.getOptionValue("samplesFile")));
				VariantAlleleFrequency variantAlleleFrequency = new VariantAlleleFrequency(output, inputCSV, samples, line.getOptionValue("region"), threads);
				if(line.hasOption("pairs")) {
					PairMapping pairs = parsePairsFile(new File(line.getOptionValue("pairs")));
					variantAlleleFrequency.setPairs(pairs);
				}
				variantAlleleFrequency.run();
			}

			if (line.hasOption("time")) {
				System.out.print("Elapsed: ");
				System.out.println((System.currentTimeMillis() - start));
			}
		} catch (ParseException exp) {
			log.error("Parsing failed", exp);
		}
	}

	private static Map<String, File> parseSamplesFile(File samplesFile) {
		BufferedReader reader = null;
		Map<String, File> retVal = new HashMap<String, File>();
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(samplesFile), Charset.forName("ASCII")));
			while (true) {
				String line = reader.readLine();
				if (line == null) break;

				if (line.length() == 0 || line.startsWith("#")) {
					continue;
				}
				String values[] = line.split(",");
				retVal.put(values[0], new File(values[1]));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		}
		return retVal;
	}
	
	private static PairMapping parsePairsFile(File pairsFile) {
		BufferedReader reader = null;
		PairMapping retVal = new PairMapping();
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(pairsFile), Charset.forName("ASCII")));
			while (true) {
				String line = reader.readLine();
				if (line == null) break;

				if (line.length() == 0 || line.startsWith("#")) {
					continue;
				}
				String values[] = line.split(",");
				retVal.add(values[0], values[1], values[2]);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		}
		return retVal;
	}
}
