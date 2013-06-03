package ca.mcgill.genome;

import java.util.HashMap;
import java.util.Map;

public class PairMapping {
	private final Map<String, String> sampleToNormal = new HashMap<String, String>();
	private final Map<String, String> sampleToTumor  = new HashMap<String, String>();
	private final Map<String, String> tumorToSample  = new HashMap<String, String>();
	private final Map<String, String> normalToSample = new HashMap<String, String>();

	public void add(String sample, String normal, String tumor) {
		if(sampleToNormal.put(sample, normal) != null) throw new RuntimeException("Sample to normal collision");
		if(sampleToTumor.put(sample, tumor) != null) throw new RuntimeException("Sample to tumor collision");
		if(tumorToSample.put(tumor, sample) != null) throw new RuntimeException("Tumor to sample collision");
		if(normalToSample.put(normal, sample) != null) throw new RuntimeException("Normal to sample collision");
	}

	public String getNormalFromSample(String sample) {
		return sampleToNormal.get(sample);
	}
	public String getTumorFromSample(String sample) {
		return sampleToTumor.get(sample);
	}
	public String getSampleFromNormal(String normal) {
		return normalToSample.get(normal);
	}
	public String getSampleFromTumor(String tumor) {
		return tumorToSample.get(tumor);
	}
	public String getNormalFromTumor(String tumor) {
		return sampleToNormal.get(tumorToSample.get(tumor));
	}
	public String getTumorFromNormal(String normal) {
		return sampleToTumor.get(normalToSample.get(normal));
	}

}
