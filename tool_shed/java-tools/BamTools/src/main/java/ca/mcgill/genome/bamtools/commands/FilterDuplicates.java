package ca.mcgill.genome.bamtools.commands;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.Charset;
import java.util.zip.GZIPOutputStream;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.filterdups.DuplicateKmerBuilder;
import ca.mcgill.genome.bamtools.filterdups.SequenceLengthException;
import ca.mcgill.genome.bamtools.parsers.FASTQFileParser;
import ca.mcgill.genome.bamtools.parsers.FASTQPEFileParser;
import ca.mcgill.genome.bamtools.parsers.Sequence;
import ca.mcgill.genome.bamtools.parsers.SequenceFileParser;

public class FilterDuplicates extends DefaultTool {
	private final DuplicateKmerBuilder duplicateKmerBuilder = new DuplicateKmerBuilder();
	private File read1 = null;
	private File read2 = null;
	private byte qualOffset = 33;
	private String prefix = null;
	private boolean onlyOutputReadNames = false;

	@Override
	public String getCmdName() {
		return "filterdups";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--read1           Fastq of read1");
		System.out.println("\t--read2           Fastq of read2 (optional)");
		System.out.println("\t-Q                Fastq Quality offset. (default: "+ qualOffset + ")");
		System.out.println("\t--prefix          Output prefix (default original filename .dup.gz)");
		System.out.println("\t--readOutput      Only output duplicate read names");
		System.out.println("\t-k,--kmer         kmer size (default: " + duplicateKmerBuilder.getKmerSize() + ")");
		System.out.println("\t-o,--offset       Read offset to get kmer (default: "	+ duplicateKmerBuilder.getOffset() + ")");
	}

	@Override
	public int run(String[] args) {
		try {
			for (int idx = 1; idx < args.length; idx++) {
				if (args[idx].equals("--read1")) {
					idx++;
					read1 = new File(args[idx]);
				} else if (args[idx].equals("--read2")) {
					idx++;
					read2 = new File(args[idx]);
				} else if (args[idx].equals("-Q")) {
					idx++;
					qualOffset = Byte.parseByte(args[idx]);
				} else if (args[idx].equals("--readOutput")) {
					onlyOutputReadNames = true;
				} else if (args[idx].equals("--prefix")) {
					idx++;
					prefix = args[idx];
				} else if (args[idx].equals("--kmer") || args[idx].equals("-k")) {
					idx++;
					duplicateKmerBuilder.setKmerSize(Integer.parseInt(args[idx]));
				} else if (args[idx].equals("--offset")
						|| args[idx].equals("-o")) {
					idx++;
					duplicateKmerBuilder.setOffset(Integer.parseInt(args[idx]));
				}
			}

			if (read1 == null) {
				printUsage("read1 not set");
				return 1;
			}

			filterDups();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		return 0;
	}

	public void filterDups() {
		BufferedWriter read1Writer = null;
		BufferedWriter read2Writer = null;
		try {
			String read1OutputFile;
			if (prefix != null) {
				read1OutputFile = prefix + ".read1.gz";
			} else {
				read1OutputFile = read1 + ".dup.read1.gz";
			}

			String read2OutputFile = null;
			SequenceFileParser fileParser;
			long totalFileSize = read1.length();

			if (read2 != null) {
				if (prefix != null) {
					read2OutputFile = prefix + ".read2.gz";
				} else {
					read2OutputFile = read2 + ".dup.read2.gz";
				}

				totalFileSize += read2OutputFile.length();
				fileParser = new FASTQPEFileParser(read1, read2, qualOffset);
				duplicateKmerBuilder.setPairedSequence(true);
			} else {
				fileParser = new FASTQFileParser(read1, qualOffset);
			}

			long lastPercent = -1;
			System.out.println("Build kmer ["
					+ duplicateKmerBuilder.getKmerSize() + "] Hash");
			while (true) {
				Sequence sequence = fileParser.nextSequence();
				if (sequence == null)
					break;
				if (sequence.getSequence().length() <= (duplicateKmerBuilder
						.getKmerSize() + duplicateKmerBuilder.getOffset())) {
					continue;
				}

				duplicateKmerBuilder.processSequence(sequence);

				long percent = fileParser.getByteCount() * 100l / totalFileSize;
				if (percent != lastPercent) {
					lastPercent = percent;
					System.out.print("\rFile read: " + percent + "%");
				}
			}
			System.out.println();
			fileParser.close();

			double nbKmers = duplicateKmerBuilder.getNbKnownSequencesFound()
					.size();
			double total = duplicateKmerBuilder.getTotalNbReads();
			System.out.println("Dup: " + (100.0 - (nbKmers * 100.0d / total))
					+ '%');

			// Second pass
			fileParser = new FASTQFileParser(read1, qualOffset);
			totalFileSize = read1.length();
			read1Writer = new BufferedWriter(
					new OutputStreamWriter(new GZIPOutputStream(
							new FileOutputStream(read1OutputFile)),
							Charset.forName("ASCII")));
			if (read2 != null) {
				FASTQFileParser filePEParser = new FASTQFileParser(read2,
						qualOffset);
				read2Writer = new BufferedWriter(new OutputStreamWriter(
						new GZIPOutputStream(new FileOutputStream(
								read2OutputFile)), Charset.forName("ASCII")));

				lastPercent = -1;
				System.out.println("Extract 'best' unique reads");
				while (true) {
					Sequence read1Sequence = fileParser.nextSequence();
					Sequence read2Sequence = filePEParser.nextSequence();
					if (read1Sequence == null)
						break;

					String sequence = read1Sequence.getSequence()
							+ read2Sequence.getSequence();
					if (sequence.length() <= (duplicateKmerBuilder
							.getKmerSize() + duplicateKmerBuilder.getOffset())) {
						if (onlyOutputReadNames) {
							read1Writer.append(read1Sequence.getHeader());
							read1Writer.newLine();
							read2Writer.append(read2Sequence.getHeader());
							read2Writer.newLine();
						}
						continue;
					}

					int qualitySum = 0;
					for (int idx = 0; idx < read1Sequence.getBaseQualities().length; idx++) {
						qualitySum += read1Sequence.getBaseQualities()[idx];
					}
					for (int idx = 0; idx < read2Sequence.getBaseQualities().length; idx++) {
						qualitySum += read2Sequence.getBaseQualities()[idx];
					}

					boolean isWritable;
					try {
						isWritable = duplicateKmerBuilder
								.getKmerCountAndMark(
										sequence,
										qualitySum
												/ (read1Sequence
														.getBaseQualities().length + read2Sequence
														.getBaseQualities().length));
					} catch (SequenceLengthException e) {
						throw new RuntimeException(
								"Kmer too long, or offset to big: "
										+ read1Sequence.getHeader());
					}
					if (isWritable && !onlyOutputReadNames) {
						read1Writer.append('@').append(
								read1Sequence.getHeader());
						read1Writer.newLine();
						read1Writer.append(read1Sequence.getSequence());
						read1Writer.newLine();
						read1Writer.append('+').append(
								read1Sequence.getHeader());
						read1Writer.newLine();
						for (int idx = 0; idx < read1Sequence
								.getBaseQualities().length; idx++) {
							read1Writer.append((char) (read1Sequence
									.getBaseQualities()[idx] + 33));
						}
						read1Writer.newLine();
						read2Writer.append('@').append(
								read2Sequence.getHeader());
						read2Writer.newLine();
						read2Writer.append(read2Sequence.getSequence());
						read2Writer.newLine();
						read2Writer.append('+').append(
								read2Sequence.getHeader());
						read2Writer.newLine();
						for (int idx = 0; idx < read2Sequence
								.getBaseQualities().length; idx++) {
							read2Writer.append((char) (read2Sequence
									.getBaseQualities()[idx] + 33));
						}
						read2Writer.newLine();
					} else if (!isWritable && onlyOutputReadNames) {
						read1Writer.append(read1Sequence.getHeader());
						read1Writer.newLine();
						read2Writer.append(read2Sequence.getHeader());
						read2Writer.newLine();
					}

					long percent = fileParser.getByteCount() * 100l
							/ totalFileSize;
					if (percent != lastPercent) {
						lastPercent = percent;
						System.out.print("\rFile read: " + percent + "%");
					}
				}
				System.out.println();
				read1Writer.close();
				read2Writer.close();
				fileParser.close();
				filePEParser.close();
			} else {
				lastPercent = -1;
				System.out.println("Extract 'best' unique reads");
				while (true) {
					Sequence read1Sequence = fileParser.nextSequence();
					if (read1Sequence == null)
						break;

					String sequence = read1Sequence.getSequence();
					if (sequence.length() <= (duplicateKmerBuilder
							.getKmerSize() + duplicateKmerBuilder.getOffset())) {
						if (onlyOutputReadNames) {
							read1Writer.append(read1Sequence.getHeader());
							read1Writer.newLine();
						}
						continue;
					}

					int qualitySum = 0;
					for (int idx = 0; idx < read1Sequence.getBaseQualities().length; idx++) {
						qualitySum += read1Sequence.getBaseQualities()[idx];
					}

					boolean isWritable;
					try {
						isWritable = duplicateKmerBuilder
								.getKmerCountAndMark(
										sequence,
										qualitySum
												/ read1Sequence
														.getBaseQualities().length);
					} catch (SequenceLengthException e) {
						throw new RuntimeException(
								"Kmer too long, or offset to big: "
										+ read1Sequence.getHeader());
					}

					if (isWritable && !onlyOutputReadNames) {
						read1Writer.append('@').append(
								read1Sequence.getHeader());
						read1Writer.newLine();
						read1Writer.append(read1Sequence.getSequence());
						read1Writer.newLine();
						read1Writer.append('+').append(
								read1Sequence.getHeader());
						read1Writer.newLine();
						for (int idx = 0; idx < read1Sequence
								.getBaseQualities().length; idx++) {
							read1Writer.append((char) (read1Sequence
									.getBaseQualities()[idx] + 33));
						}
						read1Writer.newLine();
					} else if (!isWritable && onlyOutputReadNames) {
						read1Writer.append(read1Sequence.getHeader());
						read1Writer.newLine();
					}

					long percent = fileParser.getByteCount() * 100l
							/ totalFileSize;
					if (percent != lastPercent) {
						lastPercent = percent;
						System.out.print("\rFile read: " + percent + "%");
					}
				}
				System.out.println();
				fileParser.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (read1Writer != null) {
				try {
					read1Writer.close();
				} catch (IOException e) {
					// oh well
				}
			}
			if (read2Writer != null) {
				try {
					read2Writer.close();
				} catch (IOException e) {
					// oh well
				}
			}
		}
	}

	public static void printGC() {
		runGC();
		System.out.println("Mem: " + (usedMemory() / 1024 / 1024));
	}

	private static void runGC() {
		// It helps to call Runtime.gc()
		// using several method calls:
		for (int r = 0; r < 4; ++r)
			_runGC();
	}

	private static void _runGC() {
		long usedMem1 = usedMemory(), usedMem2 = Long.MAX_VALUE;
		for (int i = 0; (usedMem1 < usedMem2) && (i < 500); ++i) {
			s_runtime.runFinalization();
			s_runtime.gc();
			Thread.yield();

			usedMem2 = usedMem1;
			usedMem1 = usedMemory();
		}
	}

	private static long usedMemory() {
		return s_runtime.totalMemory() - s_runtime.freeMemory();
	}

	private static final Runtime s_runtime = Runtime.getRuntime();
}
