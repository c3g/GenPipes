package ca.mcgill.genome.bamtools.depth;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntLongMap;
import gnu.trove.map.hash.TIntLongHashMap;

import java.io.PrintStream;
import java.text.DecimalFormat;

public class DepthInterval {
	private final long nbBasesAtDepthBins[];
	private final String name;
	private final String chromosome;
	private final int start;
	private final int end;
	private long totalBaseCoverage; // default 0
	private long totalNbAccessibleBases;  // default 0
	private long totalNbBases;  // default 0
	private long totalNbGCBases;  // default 0
	private long nbReads;  // default 0
	private long nbProperlyPaired;  // default 0

	public DepthInterval(int maxNbBasesAtDepthBins, String name) {
		this(maxNbBasesAtDepthBins, name, null, -1, -1);
	}

	public DepthInterval(int maxNbBasesAtDepthBins, String name, String chromosome, int start, int end) {
		this.nbBasesAtDepthBins = new long[maxNbBasesAtDepthBins];  // default 0
		this.name = name;
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
	}

	public void printReportHeader(TIntArrayList summaryCoverageThresholds, PrintStream out) {
		out.print("IntervalName\tStart\tStop\tNbReads\tNbProperPairs\tTotalNbCoveredBases\tTotalCoverage\tMeanCoverage\tQ75\tQ50\tQ25\tGC%");
		for(int idx=0; idx < summaryCoverageThresholds.size(); idx++) {
			int summaryCoverageThreshold = summaryCoverageThresholds.get(idx);
			out.print("\tPctBasesCoveredAt");
			out.print(summaryCoverageThreshold);
			out.print('x');
		}
	}

	public void printReport(TIntArrayList summaryCoverageThresholds, PrintStream out) {
		int q75=-1;
		int q50=-1;
		int q25=-1;
		long nbBases = 0;
		double totalNbCoveredBasesFLT = (double)totalNbAccessibleBases;
		for(int idx=0; idx < nbBasesAtDepthBins.length; idx++) {
			nbBases += nbBasesAtDepthBins[idx];
			if(q25 == -1 && nbBases >= totalNbCoveredBasesFLT*0.25) {
				q25 = idx;
			}
			if(q50 == -1 && nbBases >= totalNbCoveredBasesFLT*0.5) {
				q50 = idx;
			}
			if(q75 == -1 && nbBases >= totalNbCoveredBasesFLT*0.75) {
				q75 = idx;
				break; // Once here, it's over
			}
		}

		DecimalFormat formatter = new DecimalFormat("#0.00");
		out.print(name);
		out.print('\t');
		out.print(start);
		out.print('\t');
		out.print(end);
		out.print('\t');
		out.print(nbReads);
		out.print('\t');
		out.print(nbProperlyPaired);
		out.print('\t');
		out.print(totalNbAccessibleBases);
		out.print('\t');
		out.print(totalBaseCoverage);
		out.print('\t');
		out.print(formatter.format((double)totalBaseCoverage/totalNbCoveredBasesFLT));
		out.print('\t');
		out.print(q75);
		out.print('\t');
		out.print(q50);
		out.print('\t');
		out.print(q25);
		out.print('\t');
		if(totalNbBases > 0) {
			out.print(formatter.format((double)totalNbGCBases/(double)totalNbBases));
		}
		else {
			out.print("NA");
		}
		
		TIntLongMap coverageAt = new TIntLongHashMap();
		for(int idx=0; idx < summaryCoverageThresholds.size(); idx++) {
			int summaryCoverageThreshold = summaryCoverageThresholds.get(idx);
			coverageAt.put(summaryCoverageThreshold, 0l);
		}

		for(int idx=0; idx < nbBasesAtDepthBins.length; idx++) {
			for(int key : coverageAt.keys()) {
				if(idx >= key) {
					coverageAt.adjustValue(key, nbBasesAtDepthBins[idx]);
				}
			}
		}

		for(int key : coverageAt.keys()) {
			out.print('\t');
			out.print(formatter.format((double)coverageAt.get(key) / (double)totalNbAccessibleBases * 100.0));
		}
	}
	
	public void add(DepthInterval interval) {
		totalNbBases += interval.totalNbBases;
		totalNbGCBases += interval.totalNbGCBases;
		totalBaseCoverage += interval.totalBaseCoverage;
		totalNbAccessibleBases += interval.totalNbAccessibleBases;
		nbReads += interval.nbReads;
		nbProperlyPaired += interval.nbProperlyPaired;
		for(int idx=0; idx < interval.getNbBasesAtDepthBins().length;idx++) {
			nbBasesAtDepthBins[idx] += interval.getNbBasesAtDepthBins()[idx];
		}
	}

	public void incrementNbReads() {
		nbReads += 1l;
	}

	public void incrementNbProperlyPaired() {
		nbProperlyPaired += 1l;
	}

	public void incrementNbAccessibleBases() {
		totalNbAccessibleBases += 1l;
	}

	public void incrementNbBases() {
		totalNbBases += 1l;
	}

	public void incrementNbGCBases() {
		totalNbGCBases += 1l;
	}

	public void incrementDepthBin(int depthToIncrement) {
		totalBaseCoverage += depthToIncrement;
		if(depthToIncrement >= nbBasesAtDepthBins.length) {
			depthToIncrement = nbBasesAtDepthBins.length-1;
		}
		nbBasesAtDepthBins[depthToIncrement]++;
	}

	public long[] getNbBasesAtDepthBins() {
		return nbBasesAtDepthBins;
	}

	public String getName() {
		return name;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
}
