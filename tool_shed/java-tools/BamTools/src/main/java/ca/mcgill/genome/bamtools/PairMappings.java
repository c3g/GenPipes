package ca.mcgill.genome.bamtools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;

public class PairMappings {
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

	public static PairMappings parsePairsFile(File pairsFile) {
		BufferedReader reader = null;
		PairMappings retVal = new PairMappings();
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
