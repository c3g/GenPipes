package ca.mcgill.genome.bamtools.commands;

import gnu.trove.list.array.TIntArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.depth.ComputeDepth;
import ca.mcgill.genome.bamtools.depth.DepthInterval;

public class DepthOfCoverage extends DefaultTool {
	private File inputBAM;
	private int maxDepth = 1000;
	private boolean ommitN = false;
	private TIntArrayList summaryCoverageThresholds;
	private File refFasta;
	private byte threads = 1;
	private int minMappingQuality = 0;
	private int minBaseQuality = 0;
	private List<DepthInterval> intervals;
	private boolean computeGC = false;
	private int binsize = 0;

	public DepthOfCoverage() {
		super();
		summaryCoverageThresholds = new TIntArrayList();
		summaryCoverageThresholds.add(10);
		summaryCoverageThresholds.add(100);
		intervals = new ArrayList<DepthInterval>();
	}

	@Override
	public String getCmdName() {
		return "depthofcoverage";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--bam                        Input BAM");
		System.out.println("\t--gc                         Compute GC content. (Needs --ref)");
		System.out.println("\t--maxDepth                   Maximum depth to compute. The higher the value, the more RAM is needed. (default: " + maxDepth + ")");
		System.out.println("\t--minMappingQuality          Only use reads with a mapping quality higher than this value (default: " + minMappingQuality + ")");
		System.out.println("\t--minBaseQuality             Only count bases with a base quality higher than this value. (default: " + minBaseQuality + ")");
		System.out.println("\t--ommitN                     Don't coun't N bases. Needs reference.");
		System.out.println("\t--ref                        Indexed reference genome");
		System.out.print("\t--summaryCoverageThresholds  Compute percentage of bases covered at given values. (default: ");
		for (int idx = 0; idx < summaryCoverageThresholds.size(); idx++) {
			if (idx > 0) {
				System.out.print(',');
			}
			System.out.print(summaryCoverageThresholds.get(idx));
		}
		System.out.println(')');
		System.out.println("\t--threads                    Threads to use. (default: " + threads + ")");
		System.out.println("\t--intervals                  Intervals in bed format. (Optional)");
		System.out.println("\t--binsize                    Builds binsize intervals on each chromosome (Optional).");
	}

	@Override
	public int run(String[] args) {
		File intervalFile = null;

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--bam")) {
				idx++;
				inputBAM = new File(args[idx]);
			} else if (args[idx].equals("--gc")) {
				idx++;
				computeGC = true;
			} else if (args[idx].equals("--intervals")) {
				idx++;
				intervalFile = new File(args[idx]);
			} else if (args[idx].equals("--maxDepth")) {
				idx++;
				maxDepth = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--minMappingQuality")) {
				idx++;
				minMappingQuality = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--minBaseQuality")) {
				idx++;
				minBaseQuality = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--ommitN")) {
				ommitN = true;
			} else if (args[idx].equals("--ref")) {
				idx++;
				refFasta = new File(args[idx]);
			} else if (args[idx].equals("--threads")) {
				idx++;
				threads = Byte.parseByte(args[idx]);
			} else if (args[idx].equals("--binsize")) {
				idx++;
				binsize = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--summaryCoverageThresholds")) {
				idx++;
				String values[] = args[idx].split(",");
				summaryCoverageThresholds.clear();
				for (String value : values) {
					summaryCoverageThresholds.add(Integer.parseInt(value));
				}
			}
		}

		if (inputBAM == null) {
			printUsage("Missing inputBAM");
			return 1;
		}
		if (binsize > 0 && intervalFile != null) {
			printUsage("Only one of binsize or intervals can be used");
			return 1;
		}
		if (refFasta != null && !computeGC && !ommitN) {
			System.err.println("WARN: Turning ref off since neither gc or ommitN were passed");
			refFasta = null;
		}
		if (computeGC && refFasta == null) {
			printUsage("gc needs ref");
			return 1;
		}
		if (ommitN && refFasta == null) {
			printUsage("ommitN needs ref");
			return 1;
		}

		if (intervalFile != null) {
			intervals = parseIntervals(maxDepth, intervalFile);
		} else {
			SAMFileReader reader = new SAMFileReader(inputBAM);
			reader.setValidationStringency(ValidationStringency.SILENT);
			try {
				SAMFileHeader header = reader.getFileHeader();
				if (!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate) || !reader.hasIndex()) {
					System.err.println("BAM must be coordinate sorted and be indexed");
					return 2;
				}

				SAMSequenceDictionary dictionary = header.getSequenceDictionary();
				for (SAMSequenceRecord record : dictionary.getSequences()) {
					if (binsize > 0) {
						for (int idx = 1; idx <= record.getSequenceLength(); idx += binsize) {
							int end = idx + binsize;
							if (end > record.getSequenceLength()) {
								end = record.getSequenceLength();
							}
							intervals.add(new DepthInterval(maxDepth, record.getSequenceName() + ':' + idx + '-' + end, record.getSequenceName(), idx, end));
						}
					} else {
						intervals.add(new DepthInterval(maxDepth, record.getSequenceName(), record.getSequenceName(), 1, record.getSequenceLength()));
					}
				}
			} finally {
				reader.close();
			}
		}

		computeCoverage();
		return 0;
	}

	private List<DepthInterval> parseIntervals(int maxNbBasesAtDepthBins, File bedFile) {
		BufferedReader reader = null;
		List<DepthInterval> retVal = new ArrayList<DepthInterval>();
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(bedFile), Charset.forName("ASCII")));
			while (true) {
				String line = reader.readLine();
				if (line == null) break;

				if (line.length() == 0) {
					continue;
				}
				String values[] = line.split("\t");
				// The last +1 is hack unexplained yet
				retVal.add(new DepthInterval(maxNbBasesAtDepthBins, values[0], values[0], Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2]) + 1));
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

	public void computeCoverage() {
		ExecutorService bamCoverageExecutor = Executors.newFixedThreadPool(threads);

		try {
			List<Future<DepthInterval>> futureIntervals = new ArrayList<Future<DepthInterval>>();
			for (DepthInterval interval : intervals) {
				Future<DepthInterval> future = bamCoverageExecutor.submit(new ComputeDepth(inputBAM, refFasta, interval, ommitN, computeGC, minMappingQuality, minBaseQuality));
				futureIntervals.add(future);
			}

			DepthInterval globalInterval = new DepthInterval(maxDepth, "Total");
			for (Future<DepthInterval> future : futureIntervals) {
				DepthInterval interval = future.get();
				globalInterval.add(interval);
			} // for futures

			globalInterval.printReportHeader(summaryCoverageThresholds, System.out);
			System.out.println();
			globalInterval.printReport(summaryCoverageThresholds, System.out);
			System.out.println();
			for (DepthInterval interval : intervals) {
				interval.printReport(summaryCoverageThresholds, System.out);
				System.out.println();
			}
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} catch (ExecutionException e) {
			throw new RuntimeException(e);
		} finally {
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
