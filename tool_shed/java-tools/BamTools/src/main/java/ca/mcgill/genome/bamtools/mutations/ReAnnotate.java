package ca.mcgill.genome.bamtools.mutations;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import ca.mcgill.genome.bamtools.parsers.GeneMapper;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;

public class ReAnnotate {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		GeneMapper mapper = new GeneMapper(new File(args[0]));
		File reMap = new File(args[1]);
		File csv = new File(args[2]);
//		Config snpEffConfig = new Config(args[4], args[3]);
//		SnpEffectPredictor predictor =  snpEffConfig.loadSnpEffectPredictor();
		
		Map<String, CSVEntry> entries = parseCSV(csv, mapper);
		
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(reMap), "ASCII"), 1*1024*1024);
			String line = reader.readLine();
			System.out.print(line);
			System.out.println("\tRefSeqEntries");
			while(true) {
				line = reader.readLine();
				if(line == null)
					break;

				String values[] = line.split("\t");
				System.out.print(line);
				System.out.print('\t');
				if(!entries.containsKey(values[1]+'-'+values[2])) {
					System.out.println("N/A");
				}
				else {
					System.out.println(entries.get(values[1]+'-'+values[2]).refSeqTranscriptEffects);
				}
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// Oh well
					e.printStackTrace();
				}
			}
		}
	}

	public static Map<String, CSVEntry> parseCSV(File csv, GeneMapper mapper) {
		Map<String, CSVEntry> retVal = new HashMap<String, CSVEntry>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(csv), "ASCII"), 1*1024*1024);
			while(true) {
				String line = reader.readLine();
				if(line == null)
					break;

				String values[] = line.split("\t");
				String effectStrings[] = values[values.length-1].split(",");
				CSVEntry entry = new CSVEntry();
				entry.chr = values[0];
				entry.pos = Integer.parseInt(values[1]);
				entry.ref = values[2];
				entry.alt = values[3];
				for (String effectString : effectStrings) {
					VcfEffect effect = new VcfEffect(effectString, FormatVersion.FORMAT_SNPEFF_2);
					if(mapper.getTranscriptId(effect.getTranscriptId()) != null) {
						if(entry.refSeqTranscriptEffects.length() > 0)
							entry.refSeqTranscriptEffects += " , ";
						entry.refSeqTranscriptEffects +=  effect.getTranscriptId() + ":" + effect.getEffect() + ':' + effect.getAa();
					}
				}
				retVal.put(entry.chr + '-' + entry.pos, entry);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// Oh well
					e.printStackTrace();
				}
			}
		}
		
		return retVal;
	}
	
	public static class CSVEntry {
		public String chr;
		public int pos;
		public String ref;
		public String alt;
		public String refSeqTranscriptEffects = "";
	}
}
