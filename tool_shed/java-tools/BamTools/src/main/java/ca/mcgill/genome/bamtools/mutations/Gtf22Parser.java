package ca.mcgill.genome.bamtools.mutations;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import ca.mcgill.genome.bamtools.parsers.GeneMapper;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Interval;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * This class creates a SnpEffectPredictor from a GTF 2.2 file
 * 
 * References: http://mblab.wustl.edu/GTF22.html
 * 
 * @author pcingola
 */
public class Gtf22Parser {

	private static final String ATTRIBUTE_PATTERN_REGEX = "\\s*(\\S+)\\s+\"(.*?)\"\\s*;";
	private static final Pattern ATTRIBUTE_PATTERN = Pattern.compile(ATTRIBUTE_PATTERN_REGEX);
	private final Genome genome;
	private final GeneMapper mapper;
	private Region cds = new Region("CDS");
	private Region utr5 = new Region("5'UTR");
	private Region utr3 = new Region("3'UTR");
	private Region intron = new Region("intron");
	private Region intergenic = new Region("intergenic");
	private Map<String, List<Marker>> mergedGeneExons = new HashMap<String, List<Marker>>();

	public Gtf22Parser(GeneMapper mapper) {
		this.genome = new Genome();
		this.mapper = mapper;
	}

	public Genome getGenome() {
		return genome;
	}

	public List<Region> parse(File gtf) {
		BufferedReader reader = null;

		Map<String, List<Marker>> gene2exon = new HashMap<String, List<Marker>>();
		Map<String, List<Marker>> transcript2exon = new HashMap<String, List<Marker>>();
		try {
			if (gtf.getName().endsWith(".gz")) {
				reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gtf)), Charset.forName("ASCII")));
			} else {
				reader = new BufferedReader(new InputStreamReader(new FileInputStream(gtf), Charset.forName("ASCII")));
			}

			while (true) {
				String line = reader.readLine();
				if (line == null)
					break;
				parse(line, gene2exon, transcript2exon);
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
		
		
		for(String geneName : gene2exon.keySet()) {
			ArrayList<Marker> mergedExons = new ArrayList<Marker>();
			mergedGeneExons.put(geneName, new ArrayList<Marker>());
			List<Marker> exons = gene2exon.get(geneName);
			Collections.sort(exons);
			for(Marker exon : exons) {
				if(mergedExons.size() == 0) {
					mergedExons.add(exon.clone());
				}
				else {
					Interval lastExon = mergedExons.get(mergedExons.size()-1);
					if(exon.getStart() < lastExon.getEnd() && exon.getEnd() > lastExon.getEnd()) {
						lastExon.setEnd(exon.getEnd());
					}
					else {
						mergedExons.add(exon.clone());
					}
					
				}
			}
		}
		
		processTranscripts(transcript2exon);

		ArrayList<Region> retVal = new ArrayList<Region>();
		retVal.add(cds);
		retVal.add(utr5);
		retVal.add(utr3);
		retVal.add(intron);
		retVal.add(intergenic);
		return retVal;
	}

	public void processTranscripts(Map<String, List<Marker>> transcript2exon) {
		for(String transcript : transcript2exon.keySet()) {
			List<Marker> exons = transcript2exon.get(transcript);
			Collections.sort(exons);
			
			boolean inCDS=false;
			boolean in5prime=false;
			boolean in3prime=false;
			Marker exon = exons.get(0);
			if(exon.getStrand() > 0) {
				inCDS=false;
				in3prime=false;
				in5prime=true;
			}
			else {
				inCDS=false;
				in3prime=true;
				in5prime=false;
			}

			for(int idx=0; idx < exons.size(); idx++) {
				exon = exons.get(idx);
				Marker prevExon = null;
				Marker nextExon = null;

				if(exon.getId().equals("CDS"))
					continue;

				for(int nextIdx=idx-1; nextIdx > -1; nextIdx--) {
					prevExon = exons.get(nextIdx);
					if(!prevExon.getId().equals("exon")) {
						prevExon = null;
						continue;
					}
					else {
						break;
					}
				}

				for(int nextIdx=idx+1; nextIdx < exons.size(); nextIdx++) {
					nextExon = exons.get(nextIdx);
					if(nextExon.getId().equals("CDS")) {
						nextExon = null;
						continue;
					}
					else {
						break;
					}
				}

				if(exon.getStrand() > 0) {
					if(exon.getId().equals("exon")) {
						if(prevExon != null) {
							intron.add(new Marker((Marker)exon.getParent(), prevExon.getEnd(), exon.getStart(), exon.getStrand(), transcript));
						}

						if(nextExon == null) {
							utr3.add(exon.clone());
						}
						else if(nextExon.getId().equals("start_codon")){
							if(nextExon.getStart() <= exon.getEnd()) { // not on a splice site
								utr5.add(new Marker((Marker)exon.getParent(), exon.getStart(), nextExon.getStart(), exon.getStrand(), transcript));
								cds.add(new Marker((Marker)exon.getParent(), nextExon.getStart(), exon.getEnd(), exon.getStrand(), transcript));
							}
							else {
								utr5.add(new Marker((Marker)exon.getParent(), exon.getStart(), exon.getEnd(), exon.getStrand(), transcript));
							}
							inCDS=true;
							in5prime=false;
						}
						else if(nextExon.getId().equals("stop_codon")){
							if(nextExon.getEnd() <= exon.getEnd()) {
								cds.add(new Marker((Marker)exon.getParent(), exon.getStart(), nextExon.getEnd(), exon.getStrand(), transcript));
								utr3.add(new Marker((Marker)exon.getParent(), nextExon.getEnd(), exon.getEnd(), exon.getStrand(), transcript));
							}
							else {
								cds.add(new Marker((Marker)exon.getParent(), exon.getStart(), exon.getEnd(), exon.getStrand(), transcript));
							}
							inCDS=false;
							in3prime=true;
						}
						else if(nextExon.getId().equals("exon")){
							if(inCDS) {
								cds.add(exon.clone());
							}
							else if(in5prime) {
								utr5.add(exon.clone());
							}
							else if(in3prime) {
								utr3.add(exon.clone());
							}
							else {
								throw new RuntimeException("Unknown state");
							}
						}
						else {
							throw new RuntimeException("Unknown next: " + transcript + ' ' + nextExon.getId());
						}
					}
					else if(exon.getId().equals("start_codon")) {
						inCDS=true;
						in3prime=false;
						in5prime=false;
					}
					else if(exon.getId().equals("stop_codon")) {
						inCDS=false;
						in3prime=true;
						in5prime=false;
					}

				} // if strand
				else {
					if(exon.getId().equals("exon")) {
						if(prevExon != null) {
							intron.add(new Marker((Marker)exon.getParent(), prevExon.getEnd(), exon.getStart(), exon.getStrand(), transcript));
						}

						if(nextExon == null) {
							utr5.add(exon.clone());
						}
						else if(nextExon.getId().equals("start_codon")){
							if(nextExon.getEnd() <= exon.getEnd()) { // on a splice site.
								cds.add(new Marker((Marker)exon.getParent(), exon.getStart(), nextExon.getEnd(), exon.getStrand(), transcript));
								utr5.add(new Marker((Marker)exon.getParent(), nextExon.getEnd(), exon.getEnd(), exon.getStrand(), transcript));
							}
							else {
								cds.add(new Marker((Marker)exon.getParent(), exon.getStart(), exon.getEnd(), exon.getStrand(), transcript));
							}
							inCDS=false;
							in5prime=true;
						}
						else if(nextExon.getId().equals("stop_codon")){
							if(nextExon.getStart() <= exon.getEnd()) {
								utr3.add(new Marker((Marker)exon.getParent(), exon.getStart(), nextExon.getStart(), exon.getStrand(), transcript));
								cds.add(new Marker((Marker)exon.getParent(), nextExon.getStart(), exon.getEnd(), exon.getStrand(), transcript));
							}
							else {
								utr3.add(new Marker((Marker)exon.getParent(), exon.getStart(), exon.getEnd(), exon.getStrand(), transcript));
							}
							inCDS=true;
							in3prime=false;
						}
						else if(nextExon.getId().equals("exon")){
							if(inCDS) {
								cds.add(exon.clone());
							}
							else if(in5prime) {
								utr5.add(exon.clone());
							}
							else if(in3prime) {
								utr3.add(exon.clone());
							}
							else {
								throw new RuntimeException("Unknown state");
							}
						}
						else {
							throw new RuntimeException("Unknown next: " + transcript + ' ' + nextExon.getId());
						}
					}
					else if(exon.getId().equals("start_codon")) {
						inCDS=false;
						in3prime=false;
						in5prime=true;
					}
					else if(exon.getId().equals("stop_codon")) {
						inCDS=true;
						in3prime=false;
						in5prime=false;
					}
				}
			}
		}
	}

	public Map<String, List<Marker>> getMergedGeneExons() {
		return mergedGeneExons;
	}

	/**
	 * Read and parse GTF file
	 */
	protected boolean parse(String line, Map<String, List<Marker>> gene2exon, Map<String, List<Marker>> transcript2exon) {
		String fields[] = line.split("\t");

		// Ommit headers
		if (fields.length <= 6)
			return false;

		String type = fields[2];

		// Parse fields
		String chromo = fields[0];
		String source = fields[1];
		int start = Gpr.parseIntSafe(fields[3]) - 1;
		int end = Gpr.parseIntSafe(fields[4]) - 1;
		int strand = (fields[6].equals("-") ? -1 : +1);
		// int frame = (fields[7].equals(".") ? -1 :
		// Gpr.parseIntSafe(fields[7]));

		String geneId = "", transcriptId = "";
		String geneName = null;

		// Is it protein coding?
		// boolean proteinCoding = source.equals("protein_coding");
		String geneBioType = "";
		String trBioType = "";

		// Parse attributes
		if (fields.length >= 8) {
			HashMap<String, String> attrMap = parseAttributes(fields[8]);

			// Get gene and transcript ID
			geneId = attrMap.get("gene_id");
			transcriptId = attrMap.get("transcript_id");
			geneName = attrMap.get("gene_name");

			geneBioType = attrMap.get("gene_biotype"); // Note: This is ENSEMBL
			// specific
			if (geneBioType == null)
				geneBioType = attrMap.get("gene_type"); // Note: This is GENCODE
			// specific

			trBioType = attrMap.get("transcript_type"); // Note: This is GENCODE
			// specific
		}

		// Use 'source' as bioType (ENSEMBL uses this field)
		if ((trBioType == null) || trBioType.isEmpty())
			trBioType = source;

		// Transform null to empty
		if (geneId == null)
			geneId = "";
		if (transcriptId == null)
			transcriptId = "";

		//String id = type + "_" + chromo + "_" + (start + 1) + "_" + (end + 1); // Create
		// ID
		if (geneId.isEmpty()) {
			System.err.println("WARNING: Empty gene_id. This should never happen");
		}
		else {
			if(mapper.getTranscriptId(transcriptId) != null) {
				//They don't match, not the same ensembl version
				if(!geneName.equals(mapper.getGeneName(transcriptId))) {
					System.err.println("Gene names don't match: " + geneName + " " + mapper.getGeneName(transcriptId));
				}

				if(!transcript2exon.containsKey(transcriptId)) {
					transcript2exon.put(transcriptId, new ArrayList<Marker>());
				}
				if(!gene2exon.containsKey(geneName)) {
					gene2exon.put(geneName, new ArrayList<Marker>());
				}
				Chromosome chr = genome.getOrCreateChromosome(chromo);
				Marker marker = new Marker(chr, start, end, strand, type);
				transcript2exon.get(transcriptId).add(marker);
				gene2exon.get(geneName).add(marker);
			}
		}

		return true;
	}

	/**
	 * Parse attributes
	 * 
	 * @param attrs
	 * @return
	 */
	HashMap<String, String> parseAttributes(String attrStr) {
		HashMap<String, String> keyValues = new HashMap<String, String>();

		if (attrStr.length() > 0) {
			Matcher matcher = ATTRIBUTE_PATTERN.matcher(attrStr);
			while (matcher.find()) {
				if (matcher.groupCount() >= 2) {
					String key = matcher.group(1).toLowerCase();
					String value = matcher.group(2);
					keyValues.put(key, value);
				}
			}
		}

		return keyValues;
	}
}
