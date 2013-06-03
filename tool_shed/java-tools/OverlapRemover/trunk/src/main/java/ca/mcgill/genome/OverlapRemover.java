package ca.mcgill.genome;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.io.input.CountingInputStream;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class OverlapRemover {
	private static final Logger log = LoggerFactory.getLogger(OverlapRemover.class);

	public void cleanOverlap(File inputBAM, File output, int minOverlapSize, boolean printRejected) {
        SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		try {
			Set<String> regectedReads = new HashSet<String>();
			
			final long fileSize = inputBAM.length();
			CountingInputStream counter = new CountingInputStream(new FileInputStream(inputBAM));
			SAMFileReader samReader = new SAMFileReader(new BufferedInputStream(counter, 512 * 1024 * 1024));

			final SAMFileHeader header = samReader.getFileHeader();
			final SAMFileHeader outputHeader = header.clone();
			final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(outputHeader, true, output);

			int lastPercent = -1;
			long totalReads = 0;
			long rejectedReads = 0;
			for (SAMRecord record : samReader) {
				totalReads++;
				if (!record.getReadUnmappedFlag() && record.getProperPairFlag() && record.getReferenceIndex().equals(record.getMateReferenceIndex())) {
					if(regectedReads.remove(record.getReadName())) {
						rejectedReads++;
						continue;
					}
					
					List<AlignmentBlock> blocks = record.getAlignmentBlocks();
					AlignmentBlock firstBlock = blocks.get(0);
					if(firstBlock.getReferenceStart() < record.getMateAlignmentStart()) {
						AlignmentBlock lastBlock = blocks.get(blocks.size() - 1);
						int blockEnd = lastBlock.getReferenceStart() + lastBlock.getLength();
						if ((blockEnd - record.getMateAlignmentStart()) > minOverlapSize) {
							regectedReads.add(record.getReadName());
							rejectedReads++;
							if(printRejected)
								log.info("Rejecting: " + record.getReadName() + " " + lastBlock.getLength() + " " + blockEnd + " " + record.getMateAlignmentStart() + " " + minOverlapSize + " " + regectedReads.size());
						} else {
							writer.addAlignment(record);
						}
					} else {
						// if first isn't bad this one isn't bad
						writer.addAlignment(record);
					}
				} else {
					writer.addAlignment(record);
				}

				int percent = (int)(counter.getByteCount()*100l / fileSize);
				if(percent != lastPercent) {
					log.info("Progess: {}%", percent);
					lastPercent = percent;
				}
			}
			writer.close();
			samReader.close();
			log.info("Progess: 100%");
			log.info("Percent rejected: {}", (rejectedReads*100l/totalReads));
			log.info("Rejected        : " + rejectedReads + "/" + totalReads);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Options options = new Options();
		options.addOption("h", "help", false, "print this message");
		options.addOption("o", "output", true, "output image");
		options.addOption("i", "input", true, "input summary");
		options.addOption("s", "size", true, "min overlap-size");
		options.addOption("r", "reject", false, "print rejected reads");
		options.addOption("t", "time", false, "display elapsed time");

		CommandLineParser parser = new PosixParser();
		try {
			CommandLine line = parser.parse(options, args);
			long start = System.currentTimeMillis();
			if (line.hasOption("help") || !line.hasOption("output") || !line.hasOption("input")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp(OverlapRemover.class.getName(), options);
			} else {
				OverlapRemover run = new OverlapRemover();
				run.cleanOverlap(new File(line.getOptionValue("input")), new File(line.getOptionValue("output")), Integer.parseInt(line.getOptionValue("size")), line.hasOption('r'));

			}

			if (line.hasOption("time")) {
				System.out.print("Elapsed: ");
				System.out.println((System.currentTimeMillis() - start));
			}
		} catch (ParseException exp) {
			log.error("Parsing failed", exp);
		}
	}

}
