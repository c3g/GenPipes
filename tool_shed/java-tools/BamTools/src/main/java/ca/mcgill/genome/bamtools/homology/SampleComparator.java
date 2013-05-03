package ca.mcgill.genome.bamtools.homology;

import gnu.trove.map.hash.TObjectIntHashMap;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;

import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.commands.CompareFrequency.HomozygousMatch;
import ca.mcgill.genome.bamtools.commands.CompareFrequency.SampleDetails;
import ca.mcgill.genome.bamtools.homology.SampleComparator.ComparisonResults;

public class SampleComparator implements Callable<ComparisonResults> {
	private final SampleDetails refSample;
	private final SampleDetails sampleToCompare;
	private final List<AlleleCounts> sampleAlleleCounts;
	private final double hetThreshold;
	private final int minDepth;
	private final HomozygousMatch homMatch;

	public SampleComparator(SampleDetails refSample, SampleDetails sampleToCompare, List<AlleleCounts> sampleAlleleCounts, int minDepth, double hetThreshold, HomozygousMatch homMatch) {
		this.refSample = refSample;
		this.sampleToCompare = sampleToCompare;
		this.sampleAlleleCounts = sampleAlleleCounts;
		this.minDepth = minDepth;
		this.hetThreshold = hetThreshold;
		this.homMatch = homMatch;
	}

	public static List<AlleleCounts> parseAlleleCounts(File sampleCounts) {
		BufferedReader reader = null;
		List<AlleleCounts> positions = new ArrayList<AlleleCounts>();
		try {
			reader = new BufferedReader(new InputStreamReader(new BufferedInputStream(new FileInputStream(sampleCounts), 1024*1024), "ASCII"), 100 * 1024);
			reader.readLine(); // skip header
			String line;
			while ((line = reader.readLine()) != null) {
				String nextLine[] = line.split(",");
				String chromosome = nextLine[0];
				int position = Integer.parseInt(nextLine[1]);
				AlleleCounts counts = new AlleleCounts(chromosome, position);
				String alleles[] = nextLine[2].split(" ");
				TObjectIntHashMap<String> bamCounts = new TObjectIntHashMap<String>(alleles.length);
				for(String allele : alleles) {
					String parts[] = allele.split(":");
					int count = Integer.parseInt(parts[1]);
					bamCounts.put(parts[0], count);
				}
				
				counts.setBamCounts(bamCounts);
				positions.add(counts);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if(reader != null) try {
				reader.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}				
		}
		return positions;
	}

	@Override
	public ComparisonResults call() throws Exception {
		long nbTests=0;
		long nbPassed=0;

		List<AlleleCounts> sampleToCompareAlleleCounts = parseAlleleCounts(sampleToCompare.getAlleleCountsFile());
		if(sampleToCompareAlleleCounts.size() != sampleAlleleCounts.size()) {
			throw new RuntimeException("Both sample genotypes must be of the same size: " + sampleToCompare.getSample() + " " +sampleToCompare.getAlias());
		}

		for(int idx=0; idx < sampleAlleleCounts.size(); idx++) {
			AlleleCounts sample1Count = sampleAlleleCounts.get(idx);
			AlleleCounts sample2Count = sampleToCompareAlleleCounts.get(idx);

			if(sample1Count.getDepth() < minDepth || sample2Count.getDepth() < minDepth) {
				continue;
			}

			Set<String> kept1Bases = sample1Count.getBasesAboveFraction(hetThreshold);
			Set<String> kept2Bases = sample2Count.getBasesAboveFraction(hetThreshold);
			// There might be more than 3 alleles, with low depth samples, use all of them anyways.
			if(kept1Bases.size() == 0 || kept2Bases.size() == 0) {
				System.err.println("No bases to keep: "+sample1Count.getChromosome() + ':' + sample1Count.getPosition());
				continue;
			}

			if(homMatch == HomozygousMatch.NONE) {
				Set<String> refKept1Bases = new HashSet<String>(kept1Bases);
				Set<String> refKept2Bases = new HashSet<String>(kept2Bases);
				kept1Bases.removeAll(refKept2Bases);
				kept2Bases.removeAll(refKept1Bases);
				nbTests++;
				if(kept1Bases.size() == 0 && kept2Bases.size() == 0) {
					nbPassed++;
				}
			}
			else if(homMatch == HomozygousMatch.ONE) {
				if(kept1Bases.size() == 1 && kept2Bases.size() != 1) {
					//failed
					nbTests++;
				}
				else if(kept2Bases.size() == 1 && kept1Bases.size() != 1) {
					//failed
					nbTests++;
				}
				else if(kept2Bases.size() == 1 && kept1Bases.size() == 1) {
					nbTests++;
					Set<String> refKept1Bases = new HashSet<String>(kept1Bases);
					Set<String> refKept2Bases = new HashSet<String>(kept2Bases);
					kept1Bases.removeAll(refKept2Bases);
					kept2Bases.removeAll(refKept1Bases);
					if(kept1Bases.size() == 0 && kept2Bases.size() == 0) {
						nbPassed++;
					}
				}
			}
			else if(homMatch == HomozygousMatch.BOTH) {
				if(kept1Bases.size() == 1 && kept2Bases.size() == 1) {
					nbTests++;
					//It's a mirror if the removal on kept1 ==0, the same would happen on 2, they have the same base.
					kept1Bases.removeAll(kept2Bases);
					if(kept1Bases.size() == 0) {
						nbPassed++;
					}
				}
			}
		}
		
		return new ComparisonResults(nbTests, nbPassed, refSample, sampleToCompare);
	}

	public static class ComparisonResults {
		// Not thread safe
		private static DecimalFormat formatter = new DecimalFormat("#0.00");
		private final String sampleName;
		private final String sampleToCompareName;
		private final String sampleAlias;
		private final String sampleToCompareAlias;
		private final long nbTests;
		private final long nbPassed;

		public ComparisonResults(long nbTests, long nbPassed, SampleDetails refSample, SampleDetails sampleToCompare) {
			super();
			this.nbTests = nbTests;
			this.nbPassed = nbPassed;
			this.sampleName = refSample.getSample();
			this.sampleToCompareName = sampleToCompare.getSample();
			this.sampleAlias = refSample.getAlias();
			this.sampleToCompareAlias = sampleToCompare.getAlias();
		}

		public void printHeader(PrintStream out) {
			out.println("Sample:" + sampleName + " Alias: " + sampleAlias);
			out.println("CmpSample,CmpAlias,NbPassed,NbTests,MatchFraction,Outcome");
		}
		public void printResult(double matchFraction, PrintStream out) {
			double test = (double)nbPassed/(double)nbTests;
			
			String testResult = "PASSED";
			if(!sampleName.equals(sampleToCompareName) && test >= matchFraction) {
				testResult = "FAILED";
			}
			else if(sampleName.equals(sampleToCompareName) && test < matchFraction) {
				testResult = "FAILED";
			}
			
			out.print(sampleToCompareName);
			out.print(',');
			out.print(sampleToCompareAlias);
			out.print(',');
			out.print(nbPassed);
			out.print(',');
			out.print(nbTests);
			out.print(',');
			out.print(formatter.format(test));
			out.print(',');
			out.println(testResult);
		}
		
		
	}
}
