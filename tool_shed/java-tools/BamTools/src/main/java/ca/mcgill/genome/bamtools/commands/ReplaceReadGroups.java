package ca.mcgill.genome.bamtools.commands;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import ca.mcgill.genome.bamtools.DefaultTool;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;

public class ReplaceReadGroups extends DefaultTool {

	@Override
	public String getCmdName() {
		return "replaceRG";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--bam     BAM file");
		System.out.println("\t--output  Output BAM");
		System.out.println("\t--rgid    Read Group ID");
		System.out.println("\t--rglb    Read Group Library");
		System.out
				.println("\t--rgpl    Read Group platform (e.g. illumina, solid)");
		System.out
				.println("\t--rgpu    Read Group platform unit (eg. run barcode)");
		System.out.println("\t--rgsm    Read Group sample name");
		System.out.println("\t--rgcn    Read Group sequencing center name");
		System.out.println("\t--rgds    Read Group description");
	}

	@Override
	public int run(String[] args) {
		File bam = null;
		File outputBAM = null;
		String rgid = null;
		String rglb = null;
		String rgpl = null;
		String rgpu = null;
		String rgsm = null;
		String rgcn = null;
		String rgds = null;

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--bam")) {
				idx++;
				bam = new File(args[idx]);
			} else if (args[idx].equals("--output")) {
				idx++;
				outputBAM = new File(args[idx]);
			} else if (args[idx].equals("--rgid")) {
				idx++;
				rgid = args[idx];
			} else if (args[idx].equals("--rglb")) {
				idx++;
				rglb = args[idx];
			} else if (args[idx].equals("--rgpl")) {
				idx++;
				rgpl = args[idx];
			} else if (args[idx].equals("--rgpu")) {
				idx++;
				rgpu = args[idx];
			} else if (args[idx].equals("--rgsm")) {
				idx++;
				rgsm = args[idx];
			} else if (args[idx].equals("--rgcn")) {
				idx++;
				rgcn = args[idx];
			} else if (args[idx].equals("--rgds")) {
				idx++;
				rgds = args[idx];
			}
		}

		if (bam == null && outputBAM == null && rgid == null && rgpl == null
				&& rgsm == null && rgpu == null) {
			printUsage("bam, output and rgid not set");
			return 1;
		}
		SAMFileReader in = new SAMFileReader(bam);
		in.setValidationStringency(ValidationStringency.SILENT);

		// create the read group we'll be using
		SAMReadGroupRecord rg = new SAMReadGroupRecord(rgid);
		rg.setLibrary(rglb);
		rg.setPlatform(rgpl);
		rg.setSample(rgsm);
		rg.setPlatformUnit(rgpu);
		if (rgcn != null)
			rg.setSequencingCenter(rgcn);
		if (rgds != null)
			rg.setDescription(rgds);

		System.out.println(String.format("Created read group ID=%s PL=%s LB=%s SM=%s%n", rg.getId(), rg.getPlatform(), rg.getLibrary(), rg.getSample()));

		// create the new header and output file
		final SAMFileHeader inHeader = in.getFileHeader();
		final SAMFileHeader outHeader = inHeader.clone();
		List<SAMReadGroupRecord> samReadGroupRecords = inHeader.getReadGroups();
		boolean foundRG = false;
		List<SAMReadGroupRecord> newReadGroups = new ArrayList<SAMReadGroupRecord>(
				samReadGroupRecords.size());
		for (SAMReadGroupRecord samReadGroupRecord : samReadGroupRecords) {
			if (samReadGroupRecord.getId().equals(rg.getId())) {
				foundRG = true;
				newReadGroups.add(rg);
			} else {
				newReadGroups.add(samReadGroupRecord);
			}
		}

		if (!foundRG) {
			// throw new RuntimeException("Read group didn't exist: " + rgid);
			System.err.println(String.format("Read group didn't exist: %s", rgid));
		}

		outHeader.setReadGroups(newReadGroups);
		// if (SORT_ORDER != null) outHeader.setSortOrder(SORT_ORDER);

		final SAMFileWriter outWriter = new SAMFileWriterFactory()
				.makeSAMOrBAMWriter(outHeader,
						outHeader.getSortOrder() == inHeader.getSortOrder(),
						outputBAM);

		for (final SAMRecord read : in) {
			outWriter.addAlignment(read);
		}

		// cleanup
		in.close();
		outWriter.close();
		return 0;
	}
}
