package ca.mcgill.genome.bamtools.parsers;

import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.input.CountingInputStream;

import ca.mcgill.genome.bamtools.PairMappings;

import au.com.bytecode.opencsv.CSVReader;

public class ParseVariantTSV implements Closeable{
	private File inputTSV;
	private PairMappings pairMappings;
	private Set<String> samplesWithBAM;
	private CountingInputStream countInput;
	private CSVReader reader;
	private Map<String, Integer> sampleToColumn = new HashMap<String, Integer>();
	private Map<Integer, String> columnToSample = new HashMap<Integer, String>();
	private int chromosomeIdx = -1;
	private int positionIdx = -1;
	private int geneIdx = -1;
	private int refAlleleIdx = -1;
	private int altAllelesIdx = -1;

	public ParseVariantTSV(File inputTSV) {
		this(inputTSV, null, new HashSet<String>());
	}

	public ParseVariantTSV(File inputTSV, PairMappings pairMappings, Set<String> samplesWithBAM) {
		this.inputTSV = inputTSV;
		this.pairMappings = pairMappings;
		this.samplesWithBAM = samplesWithBAM;
	}

	public String[] readNext() throws IOException {
		return reader.readNext();
	}
	
	public String[] parseHeader() {
		String[] nextLine = null;
		try {
			countInput = new CountingInputStream(new BufferedInputStream(new FileInputStream(inputTSV), 16 * 1024 * 1024));
			reader = new CSVReader(new InputStreamReader(countInput, Charset.forName("ASCII")), '\t');

			nextLine = reader.readNext();
			if (nextLine == null) throw new RuntimeException("CSV is empty");

			for (int idx = 0; idx < nextLine.length; idx++) {
				nextLine[idx] = nextLine[idx].trim();
				if (nextLine[idx].equalsIgnoreCase("chromosome")) {
					chromosomeIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("position")) {
					positionIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("gene_name")) {
					geneIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("ref_allele")) {
					refAlleleIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("alt_alleles")) {
					altAllelesIdx = idx;
				} else {
					String sample = nextLine[idx];
					if (pairMappings != null) {
						String normal = pairMappings.getNormalFromSample(sample);
						if (normal != null) {
							if (samplesWithBAM.contains(normal)) {
								sampleToColumn.put(sample, idx);
								columnToSample.put(idx, sample);
							}
						}
						String tumor = pairMappings.getTumorFromSample(sample);
						if (tumor != null) {
							if (samplesWithBAM.contains(tumor)) {
								sampleToColumn.put(sample, idx);
								columnToSample.put(idx, sample);
							}
						}
					} else {
						if (samplesWithBAM.contains(sample)) {
							sampleToColumn.put(sample, idx);
							columnToSample.put(idx, sample);
						}
					}
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return nextLine;
	}

	@Override
	public void close() throws IOException {
		if(reader != null)
			reader.close();
		reader = null;
		countInput = null;
	}


	public long getByteCount() {
		return countInput.getByteCount();
	}

	public int getChromosomeIdx() {
		return chromosomeIdx;
	}

	public int getPositionIdx() {
		return positionIdx;
	}

	public int getGeneIdx() {
		return geneIdx;
	}

	public int getRefAlleleIdx() {
		return refAlleleIdx;
	}

	public int getAltAllelesIdx() {
		return altAllelesIdx;
	}

	public Map<String, Integer> getSampleToColumn() {
		return sampleToColumn;
	}

	public Map<Integer, String> getColumnToSample() {
		return columnToSample;
	}

}
