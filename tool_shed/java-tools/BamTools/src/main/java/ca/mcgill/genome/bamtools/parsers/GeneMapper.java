package ca.mcgill.genome.bamtools.parsers;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;

import au.com.bytecode.opencsv.CSVReader;

public class GeneMapper {
	private final File geneAnnotFile;

	private int ensemblGeneIdIdx = -1;
	private int ensemblTranscriptIdIdx = -1;
	private int geneNameIdx = -1;
	private int entrezGeneIdIdx = -1;
	private int hugoNomenclatureIdIdx = -1;
	private int hugoNomenclatureSymbolIdx = -1;
	
	private final Map<String, Map<Integer, String>> gene2Annot = new HashMap<String, Map<Integer, String>>();

	public GeneMapper(File geneAnnotFile) {
		this.geneAnnotFile = geneAnnotFile;
		parse();
	}

	public String getTranscriptId(String ensemblGeneId) {
		Map<Integer, String> annotation = gene2Annot.get(ensemblGeneId);
		if(annotation != null) {
			return annotation.get(ensemblTranscriptIdIdx);
		}
		return null;
	}

	public String getGeneName(String ensemblGeneId) {
		Map<Integer, String> annotation = gene2Annot.get(ensemblGeneId);
		if(annotation != null) {
			return annotation.get(geneNameIdx);
		}
		return null;
	}
	
	public String getEntrezGeneId(String ensemblGeneId) {
		Map<Integer, String> annotation = gene2Annot.get(ensemblGeneId);
		if(annotation != null) {
			return annotation.get(entrezGeneIdIdx);
		}
		return null;
	}

	public String getHUGONomenclatureId(String ensemblGeneId) {
		Map<Integer, String> annotation = gene2Annot.get(ensemblGeneId);
		if(annotation != null) {
			return annotation.get(hugoNomenclatureIdIdx);
		}
		return null;
	}

	public String getHUGONomenclatureSymbol(String ensemblGeneId) {
		Map<Integer, String> annotation = gene2Annot.get(ensemblGeneId);
		if(annotation != null) {
			return annotation.get(hugoNomenclatureSymbolIdx);
		}
		return null;
	}

	public void parse() {
		CSVReader reader = null;
		try {
			LineNumberReader lineReader = new LineNumberReader(new InputStreamReader(new FileInputStream(geneAnnotFile), Charset.forName("ASCII")));
			reader = new CSVReader(lineReader, '\t');
			parseHeader(reader.readNext());
			
			while(true) {
				String values[] = reader.readNext();
				if(values == null) {
					break;
				}
				
//				String ensemblGeneId = values[ensemblGeneIdIdx];
//				if(ensemblGeneId == null || ensemblGeneId.trim().length() == 0) {
//					System.err.println("Skipping field. Line: "+lineReader.getLineNumber());
//					continue;
//				}

				String ensemblTranscriptId = values[ensemblTranscriptIdIdx];
				if(ensemblTranscriptId == null || ensemblTranscriptId.trim().length() == 0) {
					System.err.println("Skipping field. Line: "+lineReader.getLineNumber());
					continue;
				}
				

				Map<Integer, String> annotations = new HashMap<Integer, String>();
				if(ensemblTranscriptIdIdx != -1) {
					annotations.put(ensemblTranscriptIdIdx, values[ensemblTranscriptIdIdx]);
				}
				if(geneNameIdx != -1) {
					annotations.put(geneNameIdx, values[geneNameIdx]);
				}
				if(entrezGeneIdIdx != -1) {
					annotations.put(entrezGeneIdIdx, values[entrezGeneIdIdx]);
				}
				if(hugoNomenclatureIdIdx != -1) {
					annotations.put(hugoNomenclatureIdIdx, values[hugoNomenclatureIdIdx]);
				}
				if(hugoNomenclatureSymbolIdx != -1) {
					annotations.put(hugoNomenclatureSymbolIdx, values[hugoNomenclatureSymbolIdx]);
				}
				
				gene2Annot.put(ensemblTranscriptId, annotations);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// Oh well;
					e.printStackTrace();
				}
			}
		}
	}
	
	private void parseHeader(String header[]) {
		for(int idx=0; idx < header.length; idx++) {
			if(header[idx].equals("Ensembl Gene ID")) {
				ensemblGeneIdIdx = idx;
			}
			else if(header[idx].equals("Ensembl Transcript ID")) {
				ensemblTranscriptIdIdx = idx;
			}
			else if(header[idx].equals("Associated Gene Name")) {
				geneNameIdx = idx;
			}
			else if(header[idx].equals("EntrezGene ID")) {
				entrezGeneIdIdx = idx;
			}
			else if(header[idx].equals("HGNC ID(s)")) {
				hugoNomenclatureIdIdx = idx;
			}
			else if(header[idx].equals("HGNC symbol")) {
				hugoNomenclatureSymbolIdx = idx;
			}
		}

		if(ensemblGeneIdIdx == -1) {
			throw new RuntimeException("Missing Ensembl Gene ID");
		}
		if(ensemblTranscriptIdIdx == -1) {
			throw new RuntimeException("Missing Ensembl Transcript ID");
		}
	}
}
