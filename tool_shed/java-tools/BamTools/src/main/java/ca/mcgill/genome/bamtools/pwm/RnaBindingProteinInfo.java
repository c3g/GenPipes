package ca.mcgill.genome.bamtools.pwm;

import java.io.IOException;

import ca.mcgill.genome.bamtools.parsers.CISBPParser;

public class RnaBindingProteinInfo {
	private final CISBPParser parser;
	private final String geneName;
	private final String geneId;
	private final String cisbpRnaId;
	private final String motifId;
	private final String studyId;
	private PositionWeightMatrix pwm = null;
	
	public RnaBindingProteinInfo(CISBPParser parser, String cisbpRnaId, String geneName, String geneId, String motifId, String studyId) {
		this.parser = parser;
		this.cisbpRnaId = cisbpRnaId;
		this.geneName = geneName;
		this.geneId = geneId;
		this.motifId = motifId;
		this.studyId = studyId;
	}

	public CISBPParser getParser() {
		return parser;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getCISBPRNAId() {
		return cisbpRnaId;
	}

	public String getGeneId() {
		return geneId;
	}

	public String getMotifId() {
		return motifId;
	}

	public String getStudyId() {
		return studyId;
	}

	public PositionWeightMatrix getPwm() {
		if(pwm == null) {
			try {
				pwm = parser.readPWM(motifId);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		return pwm;
	}

}
