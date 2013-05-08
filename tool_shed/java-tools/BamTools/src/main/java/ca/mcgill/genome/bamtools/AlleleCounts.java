package ca.mcgill.genome.bamtools;

import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.text.DecimalFormat;
import java.util.HashSet;
import java.util.Set;

public class AlleleCounts {
	private final static DecimalFormat doubleFormatter = new DecimalFormat("#0.00");

	private final String chromosome;
	private final int position;
	private TObjectIntHashMap<String> bamCounts;
	private TObjectIntHashMap<String> referenceBAMCounts;

	public AlleleCounts(String chromosome, int position) {
		super();
		this.chromosome = chromosome;
		this.position = position;
	}

	public int getRefDepth() {
		int depth = 0;
		for (String base : referenceBAMCounts.keySet()) {
			depth += referenceBAMCounts.get(base);
		}
		return depth;
	}

	public Set<String> getBasesAboveFraction(double fraction) {
		TObjectDoubleHashMap<String> alleleFractions = AlleleCounts.computeAlleleFractions(getBamCounts());
		Set<String> keptBases = new HashSet<String>();
		for (String base : alleleFractions.keySet()) {
			if (alleleFractions.get(base) > fraction) {
				keptBases.add(base);
			}
		}
		return keptBases;
	}

	public int getDepth() {
		int depth = 0;
		for (String base : bamCounts.keySet()) {
			depth += bamCounts.get(base);
		}
		return depth;
	}

	public StringBuilder format(boolean formatRef) {
		if(formatRef) {
			return formatCounts(referenceBAMCounts);
		}
		else {
			return formatCounts(bamCounts);
		}
	}

	private StringBuilder formatCounts(TObjectIntHashMap<String> counts) {
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for(String seq : counts.keySet()) {
			if(first) {
				first = false;
			}
			else {
				sb.append(' ');
			}
			sb.append(seq);
			sb.append(':');
			sb.append(counts.get(seq));
		}
		return sb;
	}

	public static StringBuilder formatFractionCounts(TObjectDoubleHashMap<String> baseCounts) {
		StringBuilder sb = new StringBuilder();
		boolean first = true;
		for(String base : baseCounts.keySet()) {
			if(first) {
				first = false;
			}
			else {
				sb.append(' ');
			}
			sb.append(base.toString());
			sb.append(':');
			sb.append(doubleFormatter.format(baseCounts.get(base)));
		}
		return sb;
	}
	
	public static  TObjectDoubleHashMap<String> computeAlleleFractions(TObjectIntHashMap<String> baseCounts) {
		TObjectDoubleHashMap<String> retVal = new TObjectDoubleHashMap<String>();
		
		double totalBaseCount = 0;
		for (String base : baseCounts.keySet()) {
			totalBaseCount += baseCounts.get(base);
		}

		for (String base : baseCounts.keySet()) {
			if(baseCounts.get(base) == 0) {
				retVal.put(base, 0);
			}
			else {
				retVal.put(base, (double) baseCounts.get(base) / totalBaseCount);
			}
		}
		return retVal;
	}

	public String getChromosome() {
		return chromosome;
	}
	public int getPosition() {
		return position;
	}

	public TObjectIntHashMap<String> getBamCounts() {
		return bamCounts;
	}

	public void setBamCounts(TObjectIntHashMap<String> bamCounts) {
		this.bamCounts = bamCounts;
	}

	public TObjectIntHashMap<String> getReferenceBAMCounts() {
		return referenceBAMCounts;
	}

	void setReferenceBAMCounts(TObjectIntHashMap<String> referenceBAMCounts) {
		this.referenceBAMCounts = referenceBAMCounts;
	}
}
