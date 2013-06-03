package ca.mcgill.genome.bamtools.parsers;

import java.io.BufferedInputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;

import org.apache.commons.io.input.CountingInputStream;

import au.com.bytecode.opencsv.CSVReader;

public class ParseSimpleTSV implements Closeable {
	private File inputTSV;
	private CountingInputStream countInput;
	private CSVReader reader;
	private int chromosomeIdx = -1;
	private int positionIdx = -1;
	private int refAlleleIdx = -1;
	private int altAllelesIdx = -1;
	private int sampleIdx = -1;

	public ParseSimpleTSV(File inputTSV) {
		this.inputTSV = inputTSV;
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
				if (nextLine[idx].equalsIgnoreCase("chr")) {
					chromosomeIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("pos")) {
					positionIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("ref")) {
					refAlleleIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("alt")) {
					altAllelesIdx = idx;
				} else if (nextLine[idx].equalsIgnoreCase("sample")) {
					sampleIdx = idx;
				}
			}

			if (chromosomeIdx == -1) { throw new RuntimeException("Missing chr column"); }
			if (positionIdx == -1) { throw new RuntimeException("Missing pos column"); }
			if (refAlleleIdx == -1) { throw new RuntimeException("Missing ref column"); }
			if (altAllelesIdx == -1) { throw new RuntimeException("Missing alt column"); }
			if (sampleIdx == -1) { throw new RuntimeException("Missing sample column"); }

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return nextLine;
	}

	@Override
	public void close() throws IOException {
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

	public int getRefAlleleIdx() {
		return refAlleleIdx;
	}

	public int getAltAllelesIdx() {
		return altAllelesIdx;
	}

	public int getSampleIdx() {
		return sampleIdx;
	}

}
