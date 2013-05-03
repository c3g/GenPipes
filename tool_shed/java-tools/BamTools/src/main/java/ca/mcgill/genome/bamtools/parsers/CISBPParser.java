package ca.mcgill.genome.bamtools.parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcgill.genome.bamtools.pwm.PositionWeightMatrix;
import ca.mcgill.genome.bamtools.pwm.RnaBindingProteinInfo;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * http://cisbp-rna.ccbr.utoronto.ca/index.php
 * CISBP-RNA Database: Catalog of Inferred RNA Binding Proteins
 * 
 * @author lletourn
 *
 */
public class CISBPParser {
	private final File rbpInformation;
	private final File rbpPWMDirectory;
	
	public CISBPParser(File rbpInformation, File rbpPWMDirectory) {
		this.rbpInformation = rbpInformation;
		this.rbpPWMDirectory = rbpPWMDirectory;

		if(!this.rbpPWMDirectory.isDirectory()) {
			throw new RuntimeException("Position weight Matrix Directory, isn't a directory");
		}
	}
	
	public Map<String, RnaBindingProteinInfo> parse() throws IOException {
		Map<String, RnaBindingProteinInfo> study2rbpsInformation = new HashMap<String,RnaBindingProteinInfo>();
		
		BufferedReader reader = null;
		try {
			reader = Gpr.reader(rbpInformation.toString());
			while(true) {
				String line = reader.readLine();
				if(line == null) {
					break;
				}
				String values[] = Gpr.split(line, '\t');
				String studyId = values[13];
				RnaBindingProteinInfo info = new RnaBindingProteinInfo(this, values[0], values[6], values[5], values[3], studyId);
				study2rbpsInformation.put(studyId, info);
			}
			return study2rbpsInformation;
		}
		finally{
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					throw new RuntimeException("Can't close reader");
				}
			}
		}
	}

	public PositionWeightMatrix readPWM(String motifId) throws IOException {
		BufferedReader reader = null;
		File pwmFile = new File(rbpPWMDirectory, motifId+".txt");
		if(!pwmFile.exists()) {
			throw new IOException("File '"+pwmFile+"' doesn't exist");
		}

		try {
			reader = Gpr.reader(pwmFile.toString());
			String line = reader.readLine();
			if(line == null) {
				return null;
			}
			String bases[] = Gpr.split(line, '\t');
			List<Character> chars = new ArrayList<Character>();
			for(int idx=1; idx < bases.length;idx++) {
				if(bases[idx].length() > 1)
					throw new RuntimeException("Base too long '"+pwmFile+"': " +bases[idx]);
				if(bases[idx].charAt(0) == 'U') {
					chars.add('T');
				} else {
					chars.add(bases[idx].charAt(0));
				}
			}

			PositionWeightMatrix pwm = new PositionWeightMatrix(chars);
			while(true) {
				line = reader.readLine();
				if(line == null) {
					break;
				}
				String values[] = Gpr.split(line, '\t');
				double weights[] = new double[bases.length];
				for(int idx=1; idx < values.length;idx++) {
					weights[idx-1] = Double.parseDouble(values[idx]);
				}
				pwm.addPositionWeights(weights);
			}
			return pwm;
		}
		catch(RuntimeException e) {
			throw new RuntimeException("Parsing: "+pwmFile, e);
		}
		finally{
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					throw new RuntimeException("Can't close reader");
				}
			}
		}
	}
}
