package ca.mcgill.genome.bamtools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CoordMath;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ReadCounter extends DefaultTool {
	private static final Logger log = LoggerFactory.getLogger(ReadCounter.class);
	private int minMapQ = 0;

	@Override
	public String getCmdName() {
		return "readCounter";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--minMapQ  Min mapping quality");
		System.out.println("\t--bam      BAM to count");
		System.out.println("\t--bed      BED intervals");
		System.out.println("\t--output   output file");
	}

	@Override
	public int run(String[] args) {
		File sampleBAM = null;
		File bed = null;
		File output = null;

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--minMapQ")) {
				idx++;
				minMapQ = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--bam")) {
				idx++;
				sampleBAM = new File(args[idx]);
			} else if (args[idx].equals("--bed")) {
				idx++;
				bed = new File(args[idx]);
			} else if (args[idx].equals("--output")) {
				idx++;
				output = new File(args[idx]);
			}
		}

		if (sampleBAM == null) {
			printUsage("bam not set");
			return 1;
		} else if (bed == null) {
			printUsage("bed not set");
			return 1;
		}
		computeCoverage(sampleBAM, bed, output);

		return 0;
	}

	private Map<String, List<Interval>> parseIntervals(File intervals) {
		HashMap<String, List<Interval>> chr2Intervals = new HashMap<String, List<Interval>>();

		LineNumberReader reader = null;
		try {
			reader = new LineNumberReader(new BufferedReader(new InputStreamReader(new FileInputStream(intervals))));
			while (true) {
				String line = reader.readLine();
				if (line == null) break;

				String values[] = line.split("\t");
				Interval i = new Interval(values[0], Integer.parseInt(values[1]), Integer.parseInt(values[2]));

				if (!chr2Intervals.containsKey(i.getChromosome())) {
					ArrayList<Interval> intervalList = new ArrayList<Interval>();
					chr2Intervals.put(i.getChromosome(), intervalList);
				}
				chr2Intervals.get(i.getChromosome()).add(i);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					log.error("Couldn't close file", e);
				}
			}
		}

		for (List<Interval> intervalList : chr2Intervals.values()) {
			Collections.sort(intervalList, new Comparator<Interval>() {

				/**
				 * {@inheritDoc}
				 */
				@Override
				public int compare(Interval o1, Interval o2) {

					return o1.getStart() - o2.getStart();
				}
			});
		}
		return chr2Intervals;
	}

	public void computeCoverage(File input, File intervals, File output) {
		Map<String, List<Interval>> chr2Intervals = parseIntervals(intervals);
		SAMFileReader samReader = new SAMFileReader(input);
		samReader.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord record : samReader) {
			if (!record.getReadUnmappedFlag() && record.getMappingQuality() >= minMapQ) {
				List<Interval> allIntervals = chr2Intervals.get(record.getReferenceName());
				if (allIntervals == null) continue;

				List<AlignmentBlock> alignmentBlocks = record.getAlignmentBlocks();
				AlignmentBlock first = alignmentBlocks.get(0);
				AlignmentBlock last = alignmentBlocks.get(alignmentBlocks.size() - 1);
				for (Interval intv : allIntervals) {
					if (first.getReferenceStart() < intv.getStop() && CoordMath.getEnd(last.getReferenceStart(), last.getLength()) > intv.getStart()) {
						intv.addRead();
					} else {
						if (intv.getStart() > CoordMath.getEnd(last.getReferenceStart(), last.getLength())) {
							break;
						}
					}
				}
			}
		}
		samReader.close();

		PrintWriter writer;
		try {
			writer = new PrintWriter(output);
		} catch (FileNotFoundException e) {
			throw new RuntimeException(e);
		}
		writer.println("Chr\tstart\tstop\tReadDepth");
		for (String chr : chr2Intervals.keySet()) {
			List<Interval> intervalList = chr2Intervals.get(chr);
			for (Interval i : intervalList) {
				writer.print(chr);
				writer.print('\t');
				writer.print(i.getStart());
				writer.print('\t');
				writer.print(i.getStop());
				writer.print('\t');
				writer.println(i.getNbReads());
			}
		}
		writer.close();
	}

	private static class Interval {
		private final String chromosome;
		private final int start;
		private final int stop;
		private int nbReads;

		public Interval(String chromosome, int start, int stop) {
			super();
			this.chromosome = chromosome;
			this.start = start;
			this.stop = stop;
			this.nbReads = 0;
		}

		public int getNbReads() {
			return nbReads;
		}

		public void addRead() {
			nbReads++;
		}

		public String getChromosome() {
			return chromosome;
		}

		public int getStart() {
			return start;
		}

		public int getStop() {
			return stop;
		}
	}
}
