package ca.mcgill.genome.bamtools.mutations;

import gnu.trove.map.hash.TObjectLongHashMap;

import java.util.HashMap;
import java.util.Map;

public class MutationRateTarget {
	public  static final String substitutionOrder[] = {"A>T", "A>C", "A>G", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"}; 

	private final Map<String, TObjectLongHashMap<String>> countPerTypePerSample = new HashMap<String, TObjectLongHashMap<String>>();
	private final Map<String, TObjectLongHashMap<String>> substitutionCountPerSample = new HashMap<String, TObjectLongHashMap<String>>();
	private final Map<String, TObjectLongHashMap<String>> substitutionCountPerTypePerSample = new HashMap<String, TObjectLongHashMap<String>>();

	public Map<String, TObjectLongHashMap<String>> getCountPerTypePerSample() {
		return countPerTypePerSample;
	}
	public Map<String, TObjectLongHashMap<String>> getSubstitutionCountPerSample() {
		return substitutionCountPerSample;
	}
	public Map<String, TObjectLongHashMap<String>> getSubstitutionCountPerTypePerSample() {
		return substitutionCountPerTypePerSample;
	}
	
}
