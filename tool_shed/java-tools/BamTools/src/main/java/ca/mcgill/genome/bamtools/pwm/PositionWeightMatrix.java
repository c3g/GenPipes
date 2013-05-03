package ca.mcgill.genome.bamtools.pwm;

import gnu.trove.map.TCharDoubleMap;
import gnu.trove.map.TCharIntMap;
import gnu.trove.map.hash.TCharIntHashMap;

import java.util.ArrayList;
import java.util.List;

public class PositionWeightMatrix {
	private final TCharIntMap bases = new TCharIntHashMap();
	private ArrayList<double[]> positionWeigths;

	public PositionWeightMatrix(List<Character> bases) {
		int idx=0;
		for(Character ch : bases) {
			this.bases.put(ch.charValue(), idx);
			idx++;
		}
		positionWeigths = new ArrayList<double[]>();
	}

	public int size() {
		return positionWeigths.size();
	}

	/**
	 * 
	 * @param weights
	 * @return position added
	 */
	public void addPositionWeights(double[] weights) {
		positionWeigths.add(weights);
	}

	public double score(String seq) {
		if(seq.length() != positionWeigths.size()) {
			throw new RuntimeException("Can't compute score, Sizes don't match");
		}
		double score = 1;
		int idx=0;
		for(double[] weights : positionWeigths) {
			score *= weights[bases.get(seq.charAt(idx))];
			idx++;
		}
		
		return score;
	}

	public double likelyhood(String seq, TCharDoubleMap atcgFrequencies) {
		if(seq.length() != positionWeigths.size()) {
			throw new RuntimeException("Can't compute score, Sizes don't match");
		}
		double score = 1;
		int idx=0;
		for(double[] weights : positionWeigths) {
			char base = seq.charAt(idx);
			score *= weights[bases.get(base)]/atcgFrequencies.get(base);
			idx++;
		}
		
		return score;
	}

	public double relativeScore(String seq, TCharDoubleMap atcgFrequencies) {
		if(seq.length() != positionWeigths.size()) {
			throw new RuntimeException("Can't compute score, Sizes don't match");
		}
		
		double minScore=0;
		double maxScore=0;

		for(double[] weights : positionWeigths) {
			double minLocalScore=Double.POSITIVE_INFINITY;
			double maxLocalScore=Double.NEGATIVE_INFINITY;
			
			for(char base : bases.keys()) {
				double score = Math.log(weights[bases.get(base)]/atcgFrequencies.get(base));
				if(score < minLocalScore)
					minLocalScore = score;
				if(score > maxLocalScore)
					maxLocalScore = score;
			}
			minScore += minLocalScore;
			maxScore += maxLocalScore;
		}

		
		double score = 1;
		int idx=0;
		for(double[] weights : positionWeigths) {
			char base = seq.charAt(idx);
			score += Math.log(weights[bases.get(base)]/atcgFrequencies.get(base));
			idx++;
		}
		
		return (score-minScore)/(maxScore-minScore);
	}
}
