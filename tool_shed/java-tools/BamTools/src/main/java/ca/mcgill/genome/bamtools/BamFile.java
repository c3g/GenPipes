package ca.mcgill.genome.bamtools;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class BamFile {

	final private String strand;
	final private SAMFileReader fileReader;
	final private String order;
	final private String fileName;
	
	public BamFile(String strand, String fileReader, String order) {
		this.strand = strand;
		this.fileReader = new SAMFileReader(new File(fileReader));
		this.fileReader.setValidationStringency(ValidationStringency.SILENT);
		this.order = order;
		this.fileName = fileReader;
		// TODO Auto-generated constructor stub
	}
	
	public String getStrand() {
		return this.strand;
	}
	
	public SAMFileReader getFileReader() {
		return this.fileReader;
	}

	public String getOrder() {
		return this.order;
	}
	public String getFileName() {
		return this.fileName;
	}
}
