package ca.mcgill.genome.bamtools.commands;

import gnu.trove.impl.Constants;
import gnu.trove.map.hash.TObjectIntHashMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.SampleAlleleCounts.ReadBases;
import ca.mcgill.genome.bamtools.homology.AlleleFrequencyComputer;

public class BaseFrequency extends DefaultTool {
	private boolean printElapsed = false;
	//private IndexedFastaSequenceFile refSequence = null;
	private byte threads = 1;
	private int minMappingQuality = 10;
	private int minBaseQuality = 10;

	@Override
	public String getCmdName() {
		return "basefreq";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--pos                  Input position in TSV. chr,1-basedPosition");
		System.out.println("\t--bam                  BAM file");
		System.out.println("\t--threads              Threads to use. (default: " + threads + ")");
		System.out.println("\t--printElapsed         Prints elapsed percentage");
		//System.out.println("\t--ref           Reference to add refAllele");
		System.out.println("\t--noIndex              Don't use the index");
		System.out.println("\t--minMappingQuality    Only use reads with a mapping quality higher than this value (default: " + minMappingQuality + ")");
		System.out.println("\t--minBaseQuality       Only count bases with a base quality higher than this value. (default: " + minBaseQuality + ")");
	}

	@Override
	public int run(String[] args) {
		File bamFile = null;
		File positions = null;
		boolean noIndex = false;

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--pos")) {
				idx++;
				positions = new File(args[idx]);
			} else if (args[idx].equals("--threads")) {
				idx++;
				threads = Byte.parseByte(args[idx]);
			} else if (args[idx].equals("--bam")) {
				idx++;
				bamFile = new File(args[idx]);
//			} else if (args[idx].equals("--ref")) {
//				idx++;
//				refSequence = new IndexedFastaSequenceFile(new File(args[idx]));
			} else if (args[idx].equals("--printElapsed")) {
				printElapsed = true;
			} else if (args[idx].equals("--noIndex")) {
				noIndex = true;
			} else if (args[idx].equals("--minMappingQuality")) {
				idx++;
				minMappingQuality = Integer.parseInt(args[idx]);
			} else if (args[idx].equals("--minBaseQuality")) {
				idx++;
				minBaseQuality = Integer.parseInt(args[idx]);
			}
		}

		if (bamFile == null) {
			printUsage("bam not set");
			return 1;
		}

		if (positions == null) {
			printUsage("pos not set");
			return 1;
		}

		List<AlleleCounts> positionsToTest = parsePositions(positions);
		if(noIndex)
			computeFreqWithoutIndex(positionsToTest, bamFile);
		else
			computeFreq(positionsToTest, bamFile);

		printCounts(positionsToTest);
		return 0;
	}

	public List<AlleleCounts> parsePositions(File positionsFile) {
		BufferedReader reader = null;
		Map<String, AlleleCounts> positions = new HashMap<String, AlleleCounts>();
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(positionsFile), "ASCII"), 100 * 1024);
			String line;
			while ((line = reader.readLine()) != null) {
				if(line.length() > 0 && line.charAt(0) == '#')
					continue;
				String nextLine[] = line.split("\t");
				String chromosome = nextLine[0];
				String key = chromosome+'#'+nextLine[1];
				int position = Integer.parseInt(nextLine[1]);
				positions.put(key, new AlleleCounts(chromosome, position));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			if(reader != null) try {
				reader.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}				
		}
		return new ArrayList<AlleleCounts>(positions.values());
	}

	public void computeFreqWithoutIndex(List<AlleleCounts> positionsToTest, File bamFile) {
		for(AlleleCounts alleleCounts : positionsToTest) {
			alleleCounts.setBamCounts(new TObjectIntHashMap<String>(Constants.DEFAULT_CAPACITY, Constants.DEFAULT_LOAD_FACTOR, 0));
			for (ReadBases bases : ReadBases.values()) {
				alleleCounts.getBamCounts().put(bases.toString(), 0);
			}
		}

		SAMFileReader samReader = new SAMFileReader(bamFile);
		samReader.setValidationStringency(ValidationStringency.SILENT);
		try {
			SAMFileHeader header = samReader.getFileHeader();
			if(!header.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
				throw new RuntimeException("BAM must be coordinate sorted");
			}
			final SAMSequenceDictionary dictionary = header.getSequenceDictionary();
			Collections.sort(positionsToTest, new Comparator<AlleleCounts>() {
				@Override
				public int compare(AlleleCounts o1, AlleleCounts o2) {
					int idx1 = dictionary.getSequenceIndex(o1.getChromosome());
					int idx2 = dictionary.getSequenceIndex(o2.getChromosome());
					int delta = idx1-idx2;
					if(delta !=0)
						return delta;
					
					return (int)(o1.getPosition()-o2.getPosition());
				}
			});

			int currentIdx = 0;
			for(SAMRecord record : samReader) {
				if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag() || record.getMappingQuality() < minMappingQuality)
					continue;

				int start = record.getAlignmentStart();
				int end = record.getAlignmentEnd();

				boolean skipRecord = false;
				while(currentIdx < positionsToTest.size()) {
					if(!record.getReferenceName().equals(positionsToTest.get(currentIdx).getChromosome())) {
						if(record.getReferenceIndex().intValue() < dictionary.getSequenceIndex(positionsToTest.get(currentIdx).getChromosome())) {
							skipRecord = true;
							break;
						}
						else {
							
						}
					}
					else if(positionsToTest.get(currentIdx).getPosition() >= start) {
						break;
					}
					currentIdx++;
				}
				if(currentIdx >= positionsToTest.size())
					break;
				if(skipRecord) 
					continue;
				
				if(positionsToTest.get(currentIdx).getPosition() >= start && positionsToTest.get(currentIdx).getPosition() <= end) {
					for(int tmpIdx = currentIdx; tmpIdx <  positionsToTest.size() && record.getReferenceName().equals(positionsToTest.get(tmpIdx).getChromosome()); tmpIdx++) {
						if(positionsToTest.get(tmpIdx).getPosition() <= end) {
							for (AlignmentBlock block : record.getAlignmentBlocks()) {
								if (block.getReferenceStart() <= positionsToTest.get(tmpIdx).getPosition() && positionsToTest.get(tmpIdx).getPosition() <= (block.getReferenceStart() + block.getLength())) {
									int refOffset = positionsToTest.get(tmpIdx).getPosition() - block.getReferenceStart();
									int readBaseOffset = block.getReadStart() - 1 + refOffset;
									if (record.getBaseQualities()[readBaseOffset] >= minBaseQuality) {
										byte bases[] = record.getReadBases();
										String baseStr = String.valueOf((char) bases[readBaseOffset]);
										positionsToTest.get(tmpIdx).getBamCounts().put(baseStr, (positionsToTest.get(tmpIdx).getBamCounts().get(baseStr) + 1));
									}
									break;
								}
							}
						}
						else {
							break;
						}
					}
				}
//				else {
//					if(positionsToTest.get(currentIdx).getPosition() > end) {
//						currentIdx++;
//					}
//				}
			}
		}
		finally {
			samReader.close();
		}
	}

	public void computeFreq(List<AlleleCounts> positionsToTest, File bamFile) {
		Queue<AlleleCounts> positions = new ConcurrentLinkedQueue<AlleleCounts>(positionsToTest);

		try {
			AlleleFrequencyComputer alleleFrequencyComputers[] = new AlleleFrequencyComputer[threads];
			Thread workers[] = new Thread[threads];
			ThreadGroup homologyGrp = new ThreadGroup("AlleleFrequencyComputer");
			for(int idx=0; idx < threads; idx++) {
				alleleFrequencyComputers[idx] = new  AlleleFrequencyComputer(positions, bamFile, minMappingQuality, minBaseQuality);
				workers[idx] = new Thread(homologyGrp, alleleFrequencyComputers[idx]);
				workers[idx].setDaemon(true);
				workers[idx].start();
			}

			long percent = 0;
			int wait = 1000;
			for(int idx=0; idx < workers.length; idx++) {
				Thread worker = workers[idx];
				while(worker.isAlive()) {
					Thread.sleep(wait);
					if(printElapsed) {
						long currentPercent = (positionsToTest.size()-positions.size())*100l/positionsToTest.size();
						if(percent != currentPercent) {
							percent = currentPercent;
							System.out.print("\rCompletion: " + percent + '%');
						}
					}
				}
				try {
					alleleFrequencyComputers[idx].close();
				} catch (IOException e) {
					System.err.println("Problems closing Frequency Computer:" + e);
				}
			}
			if(printElapsed) {
				System.out.println("\rCompletion: 100%");
			}
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
	}
	
	public void printCounts(List<AlleleCounts> positionsToTest) {
		System.out.println("Chr,Position,AlleleCounts,AlleleFrequency");
		for(AlleleCounts alleleCounts : positionsToTest) {
			System.out.print(alleleCounts.getChromosome());
			System.out.print(',');
			System.out.print(alleleCounts.getPosition());
			System.out.print(',');
			System.out.print(alleleCounts.format(false));
			System.out.print(',');
			System.out.println(AlleleCounts.formatFractionCounts(AlleleCounts.computeAlleleFractions(alleleCounts.getBamCounts())));
		}
	}
}
