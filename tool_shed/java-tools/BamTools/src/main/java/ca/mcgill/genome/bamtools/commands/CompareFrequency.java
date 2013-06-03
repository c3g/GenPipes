package ca.mcgill.genome.bamtools.commands;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.homology.SampleComparator;
import ca.mcgill.genome.bamtools.homology.SampleComparator.ComparisonResults;

public class CompareFrequency extends DefaultTool {
	public  enum HomozygousMatch {NONE, ONE, BOTH};

	private HomozygousMatch homMatch = HomozygousMatch.NONE;
	private int threads = 1;
	private int minDepth = 3;
	private double hetThreshold = 0.3;  
	private double matchFraction = 0.95;

	@Override
	public String getCmdName() {
		return "cmpfreq";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--het           At which threshold do we call a valid het. (default: " + hetThreshold + ")");
		System.out.println("\t--match         At what percent do we declare a match. (default: " + matchFraction + ")");
		System.out.println("\t--sample        Sample name");
		System.out.println("\t--sampleCounts  Sample counts file");
		System.out.println("\t--compareList   List of files to compare to. Format <Name,Alias,Path>. The name will be used to see if it corresponds");
		System.out.println("\t--minDepth      Minimum Depth to consider. (default: " + minDepth + ")");
		System.out.println("\t--homMatch      Use only homozygous calls. (default: " + homMatch + ")");
		System.out.println("\t                "+HomozygousMatch.NONE+" : Use Homozygous and Heterozygous calls");
		System.out.println("\t                "+HomozygousMatch.ONE +" : One of the 2 samples must be homozygous");
		System.out.println("\t                "+HomozygousMatch.BOTH+" : Both need to be homozygous");
		System.out.println("\t--threads       Threads to use. (default: "+threads+")");
	}

	@Override
	public int run(String[] args) {
		String sampleName = null;
		File sampleCountsFile = null;
		File compareListFile = null;

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--sample")) {
				idx++;
				sampleName = args[idx];
			} else if (args[idx].equals("--het")) {
				idx++;
				hetThreshold = Double.parseDouble(args[idx]);
			} else if (args[idx].equals("--match")) {
				idx++;
				matchFraction = Double.parseDouble(args[idx]);
			} else if (args[idx].equals("--sampleCounts")) {
				idx++;
				sampleCountsFile = new File(args[idx]);
			} else if (args[idx].equals("--compareList")) {
				idx++;
				compareListFile = new File(args[idx]);
			} else if (args[idx].equals("--minDepth")) {
				idx++;
				minDepth = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--threads")) {
				idx++;
				threads = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--homMatch")) {
				idx++;
				homMatch = Enum.valueOf(HomozygousMatch.class, args[idx]);
			}
		}

		if(sampleName == null || sampleCountsFile == null){
			printUsage("Both sample and sampleCounts need to be specified");
			return 1;
		}
		if(compareListFile == null){
			printUsage("Both compareList needs to be specified");
			return 1;
		}

		List<SampleDetails> compareList = parseSamplesFile(compareListFile);
		SampleDetails sample = new SampleDetails(sampleName,"", sampleCountsFile);
		compareCounts(sample, compareList);
		return 0;
	}


	public List<SampleDetails> parseSamplesFile(File compareListFile) {
		BufferedReader reader = null;
		List<SampleDetails> retVal = new ArrayList<CompareFrequency.SampleDetails>();
		try {
			reader = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(compareListFile), 1024*1024), Charset.forName("ASCII")));
			while (true) {
				String line = reader.readLine();
				if (line == null) break;

				if (line.length() == 0 || line.startsWith("#")) {
					continue;
				}
				String values[] = line.split(",");

				SampleDetails sd = new SampleDetails(values[0], values[1], new File(values[2]));
				retVal.add(sd);
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

	public void compareCounts(SampleDetails sample, List<SampleDetails> samplesToCompare) {
		ExecutorService sampleComparatorExecutor = Executors.newFixedThreadPool(threads);

		try {
			List<Future<ComparisonResults>> futures = new ArrayList<Future<ComparisonResults>>();
			List<AlleleCounts> sampleAlleleCounts = SampleComparator.parseAlleleCounts(sample.getAlleleCountsFile());
			for(SampleDetails sampleToCompare : samplesToCompare) {
				futures.add(sampleComparatorExecutor.submit(new SampleComparator(sample, sampleToCompare, sampleAlleleCounts, minDepth, hetThreshold, homMatch)));
			}
			
			boolean first = true;
			for(Future<ComparisonResults> future : futures) {
				ComparisonResults result = future.get();
				if(first) {
					first = false;
					result.printHeader(System.out);
				}
				result.printResult(matchFraction, System.out);
			}
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		} catch (ExecutionException e) {
			throw new RuntimeException(e);
		} finally {
			sampleComparatorExecutor.shutdownNow();
		}
	}
	
	public static class SampleDetails {
		private final String sample;
		private final String alias;
		private final File alleleCountsFile;

		public SampleDetails(String sample, String alias, File alleleCountsFile) {
			this.sample = sample;
			this.alias = alias;
			this.alleleCountsFile = alleleCountsFile;
		}

		public String getSample() {
			return sample;
		}

		public String getAlias() {
			return alias;
		}

		public File getAlleleCountsFile() {
			return alleleCountsFile;
		}
	}
}
