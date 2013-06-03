package ca.mcgill.genome.bamtools.commands;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.PairMappings;
import ca.mcgill.genome.bamtools.SampleAlleleCounts;
import ca.mcgill.genome.bamtools.parsers.ParseVariantTSV;

public class VariantAlleleFrequency extends DefaultTool {
	private static final Logger log = LoggerFactory.getLogger(VariantAlleleFrequency.class);
	private File output;
	private File inputCSV;
	private Map<String, SampleAlleleCounts> samplesBAM;
	private String chromosome;
	private byte threads;
	private PairMappings pairs = null;

	@Override
	public String getCmdName() {
		return "varfreq";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--output       Output file");
		System.out.println("\t--csv          Input CSV");
		System.out.println("\t--samplesFile  Samples BAMs <needed with vtsv, tsv>");
		System.out.println("\t--threads      Threads to use. (default: "+threads+")");
		System.out.println("\t--chr          Chromosome");
	}

	@Override
	public int run(String[] args) {
		File samplesToFile = null;
		samplesBAM = new HashMap<String, SampleAlleleCounts>();

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--output")) {
				idx++;
				output = new File(args[idx]);
			} else if (args[idx].equals("--csv")) {
				idx++;
				inputCSV = new File(args[idx]);
			} else if (args[idx].equals("--threads")) {
				idx++;
				threads = Byte.parseByte(args[idx]);
			} else if (args[idx].equals("--bams")) {
				idx++;
				samplesToFile = new File(args[idx]);
			} else if (args[idx].equals("--chr")) {
				idx++;
				chromosome = args[idx];
			}
		}

		if(samplesToFile == null) {
			printUsage("Missing SamplesToFile not set");
			return 1;
		}

		Map<String, File> samplesFile = parseSamplesFile(samplesToFile);
		for(String sampleName : samplesFile.keySet()) {
			samplesBAM.put(sampleName, new SampleAlleleCounts(null, samplesFile.get(sampleName), 15, 0));
		}
		computeFreq();
		return 0;
	}

	public void setPairs(PairMappings pairs) {
		this.pairs = pairs;
	}

	public void computeFreq() {
		ExecutorService bamCoverageExecutor = Executors.newFixedThreadPool(threads);

		PrintWriter writer = null;
		ParseVariantTSV parseVariantTSV = null;
		try {
			long inputCSVSize = inputCSV.length();
			writer = new PrintWriter(output, "ASCII");
			parseVariantTSV = new ParseVariantTSV(inputCSV, pairs, samplesBAM.keySet());
			String nextLine[] = parseVariantTSV.parseHeader();

			for(int idx=0; idx < nextLine.length; idx++) {
				if(idx > 0) {
					writer.print('\t');
				}
				String sample = parseVariantTSV.getColumnToSample().get(idx);
				if(sample != null) {
					writer.print(sample + "allele counts");
					writer.print('\t');
				}
				writer.print(nextLine[idx]);
			}
			writer.println();
			
			Map<String, Future<AlleleCounts>> sampleFutures = new HashMap<String, Future<AlleleCounts>>(parseVariantTSV.getSampleToColumn().size());
			long prevPercent=0;
			while ((nextLine = parseVariantTSV.readNext()) != null) {
				long percent = parseVariantTSV.getByteCount()*100l/inputCSVSize;
				if(percent != prevPercent) {
					log.info("Completed: {}%", percent);
					prevPercent = percent;
				}
				sampleFutures.clear();
				String chromosome = nextLine[parseVariantTSV.getChromosomeIdx()];
				// Skip chromosomes that don't match if asked for.
				if(this.chromosome != null && !this.chromosome.equals(chromosome))
					continue;

				int position = Integer.parseInt(nextLine[parseVariantTSV.getPositionIdx()]);

				for(String sampleName : parseVariantTSV.getSampleToColumn().keySet()) {
					if(pairs != null) {
						String normal = pairs.getNormalFromSample(sampleName);
						Future<AlleleCounts> future;
						future = bamCoverageExecutor.submit(samplesBAM.get(normal).getCallableAlleleCountsComputer(chromosome, position));
						sampleFutures.put(normal, future);
						String tumor = pairs.getTumorFromSample(sampleName);
						future = bamCoverageExecutor.submit(samplesBAM.get(tumor).getCallableAlleleCountsComputer(chromosome, position));
						sampleFutures.put(tumor, future);
					}
					else {
						Future<AlleleCounts> future = bamCoverageExecutor.submit(samplesBAM.get(sampleName).getCallableAlleleCountsComputer(chromosome, position));
						sampleFutures.put(sampleName, future);
					}
				}
				
				Map<String, StringBuilder> toPrint = new HashMap<String, StringBuilder>();
				for(String sampleName : sampleFutures.keySet()) {
					Future<AlleleCounts> future = sampleFutures.get(sampleName);
					AlleleCounts baseCounts = future.get();
					
					if(pairs != null) {
						if(pairs.getSampleFromNormal(sampleName) != null) {
							String sample = pairs.getSampleFromNormal(sampleName);
							if(!toPrint.containsKey(sample))
								toPrint.put(sample, new StringBuilder());

							StringBuilder sb = baseCounts.format(false);
							sb.append(" , ");
							toPrint.get(sample).insert(0, sb.toString());
						}
						else if(pairs.getSampleFromTumor(sampleName) != null) {
							String sample = pairs.getSampleFromTumor(sampleName);
							if(!toPrint.containsKey(sample))
								toPrint.put(sample, new StringBuilder());

							StringBuilder sb = baseCounts.format(false);
							toPrint.get(sample).append(sb);
						}
						else {
							throw new RuntimeException("Unknown sample: " + sampleName);
						}
					}
					else {
						StringBuilder sb = new StringBuilder();
						boolean first = true;
						for(String base : baseCounts.getBamCounts().keySet()) {
							if(first) {
								first = false;
							}
							else {
								sb.append(' ');
							}
							sb.append(base.toString());
							sb.append(':');
							sb.append(baseCounts.getBamCounts().get(base));
						}
						toPrint.put(sampleName, sb);
					}
					
				} // for futures
				
				for(int idx=0; idx < nextLine.length; idx++) {
					if(idx > 0) {
						writer.print('\t');
					}
					String sample = parseVariantTSV.getColumnToSample().get(idx);
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
			if (parseVariantTSV != null) {
				try {
					parseVariantTSV.close();
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

	public static Map<String, File> parseSamplesFile(File samplesFile) {
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
	
	
}
