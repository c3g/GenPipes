package ca.mcgill.genome.bamtools.commands;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import ca.mcgill.genome.bamtools.Allele;
import ca.mcgill.genome.bamtools.BamFile;
import ca.mcgill.genome.bamtools.CGSite;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.HetSite;




public class AlleleSpecificMethyl extends DefaultTool {
	
	private byte threads;
	private HashMap<Integer,Integer> tmpHetHash = new HashMap<Integer,Integer>();

	public AlleleSpecificMethyl() {
		// TODO Auto-generated constructor stub
	}

	@Override
	public String getCmdName() {
		return "asm";
	}
	
	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}
	
	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--pos           input position file in TSV");
		System.out.println("\t--fbam1          forward BAM file phase 1");
		System.out.println("\t--rbam1          reverse BAM file phase 1");
		System.out.println("\t--fbam2          forward BAM file phase 2");
		System.out.println("\t--rbam2          reverse BAM file phase 2");
		System.out.println("\t--ref           fasta reference file");
		System.out.println("\t--threads       number of threads to use. (default: " + threads + ")");
	}

	@Override
	public int run(String[] args) {
		
		String forwardBamFile1 = null;
		String reverseBamFile1 = null;
		String forwardBamFile2 = null;
		String reverseBamFile2 = null;
		File positions = null;
		File referenceFasta = null;
		

		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--pos")) {
				idx++;
				positions = new File(args[idx]);
			} else if (args[idx].equals("--threads")) {
				idx++;
				threads = Byte.parseByte(args[idx]);
			} else if (args[idx].equals("--fbam1")) {
				idx++;
				forwardBamFile1 = args[idx].toString();
			} else if (args[idx].equals("--rbam1")) {
				idx++;
				reverseBamFile1 = args[idx].toString();
			} else if (args[idx].equals("--fbam2")) {
				idx++;
				forwardBamFile2 = args[idx].toString();
			} else if (args[idx].equals("--rbam2")) {
				idx++;
				reverseBamFile2 = args[idx].toString();
			} else if (args[idx].equals("--ref")) {
				idx++;
				referenceFasta = new File(args[idx]);
			}
		}

		if (forwardBamFile1 == null || reverseBamFile1 == null || forwardBamFile2 == null || reverseBamFile2 == null || referenceFasta == null) {
			printUsage("bam files not set");
			return 1;
		} else if (positions == null) {
			printUsage("pos not set");
			return 1;
		}

		BamFile forwardBamFile1SAMReader = new BamFile("forward", forwardBamFile1, "first");
		BamFile forwardBamFile2SAMReader = new BamFile("forward", forwardBamFile2, "second");
		BamFile reverseBamFile1SAMReader = new BamFile("reverse", reverseBamFile1, "first");
		BamFile reverseBamFile2SAMReader = new BamFile("reverse", reverseBamFile2, "second");
		
		List<BamFile> tmpList = Arrays.asList(forwardBamFile1SAMReader,forwardBamFile2SAMReader,reverseBamFile1SAMReader,reverseBamFile2SAMReader);

		for(BamFile sfr : tmpList) {		
			if (!sfr.getFileReader().hasIndex()) { throw new RuntimeException("The bam file " + sfr + " does not have an index"); }
		}
		
		List<HetSite> positionsToTest = parsePositions(positions);
		IndexedFastaSequenceFile referenceFastaFile = null;
		try {
			referenceFastaFile = new IndexedFastaSequenceFile(referenceFasta);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		getAMR(positionsToTest, forwardBamFile1SAMReader, forwardBamFile2SAMReader, reverseBamFile1SAMReader, reverseBamFile2SAMReader, referenceFastaFile);
		return 0;
	}
	
	public List<HetSite> parsePositions(File positionsFile) {
		BufferedReader reader = null;
		List<HetSite> positions = new ArrayList<HetSite>();
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(positionsFile), "ASCII"), 100 * 1024);
			String line;
			while ((line = reader.readLine()) != null) {
				String nextLine[] = line.split("\t");
				String chromosome = nextLine[0];
				int position = new Integer(nextLine[1]);
				char allele1 = nextLine[2].charAt(0);
				char allele2 = nextLine[3].charAt(0);
				positions.add(new HetSite(chromosome, position, allele1, allele2));
				tmpHetHash.put(position, 1);
				tmpHetHash.put(position-1, 1);
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
		return positions;
	}
	
	public void getAMR(List<HetSite> positionsToTest, BamFile f1, BamFile f2, BamFile r1, BamFile r2, IndexedFastaSequenceFile referenceFastaFile) {
			
		
		//Queue<HetSite> positions = new ConcurrentLinkedQueue<HetSite>(positionsToTest);
		
		final int readLength = 100;
		HashMap<String,CGSite> globalCGSites = new HashMap<String,CGSite>();
		
		for(HetSite hetsite : positionsToTest) {
			
			String chromosome = hetsite.getChromosome();
		
			int position = hetsite.getPosition();
			//System.out.println(position);

			Character allele1 = hetsite.getAllele1();
			Character allele2 = hetsite.getAllele2();
			
			HashMap<String,CGSite> localCGSites = new HashMap<String,CGSite>();
			
			// Get CG positions from reference in readLength bp window around position

			ReferenceSequence sequence = referenceFastaFile.getSubsequenceAt(chromosome, position-readLength+1, position+readLength-1);
			String bases = null;
			try {
				bases = new String(sequence.getBases(), "ASCII");
			} catch (UnsupportedEncodingException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			//System.out.println(bases);			
			for (int i = 0; i < bases.length()-1;i++) {
				Character a = bases.charAt(i);
				Character b = bases.charAt(i+1);
				if(a.equals('C') && b.equals('G')) {
					int cgPosition = position-readLength+i+1;
				
					if(tmpHetHash.containsKey(cgPosition)) continue;
					
					String hashKey = chromosome + "," + cgPosition;
					if(!globalCGSites.containsKey(hashKey)) {
						CGSite tmpSite = new CGSite(chromosome, position);
						tmpSite.allele1.setBase(allele1);
						tmpSite.allele2.setBase(allele2);
						localCGSites.put(hashKey, tmpSite);
						globalCGSites.put(hashKey, tmpSite);
					}
				}
			}
		   
		
			List<BamFile> SAMReaderList = null;
			SAMReaderList = Arrays.asList(f1,r1,f2,r2);
			//SAMReaderList = Arrays.asList(f1);
			// Get reads overlapping het site position
			try {
				for(BamFile sfr : SAMReaderList) {
					SAMRecordIterator samIterator = null;
	
					samIterator = sfr.getFileReader().queryOverlapping(chromosome, position-1, position-1);
					//samIterator = sfr.getFileReader().queryOverlapping(chromosome, position, position);
					String bamOrder = sfr.getOrder();
					String bamStrand = sfr.getStrand();
					
					// For each read found overlapping het site
					
					while (samIterator.hasNext()) {
						
						SAMRecord record = samIterator.next();
						
						if (record.getReadUnmappedFlag() || record.getDuplicateReadFlag())
							continue;
						
						String readName = record.getReadName();
						//String readString = record.getReadString();
					
						//System.out.println(readName);
						
						//System.out.println(readString);
						
						// Get alignment blocks for each read
						
						for (AlignmentBlock block : record.getAlignmentBlocks()) {
							//System.out.println(readName);
							// Look if the block overlaps a CG site
							//System.out.println("Block ref start");
							int blockRefStart = block.getReferenceStart();
							//System.out.println(blockRefStart);
							int blockReadStart = block.getReadStart();
							//System.out.println("Block read start");
							//System.out.println(blockReadStart);
							int blockLength = block.getLength();
							//System.out.println("Block length");
							//System.out.println(blockLength);
							int refOffsetHetSite = position - blockRefStart;
							//System.out.println("refOffsetHetSite");
							//System.out.println(refOffsetHetSite);
							int readBaseOffsetHetSite = blockReadStart -1 + refOffsetHetSite;
							//System.out.println("readBaseOffsetHetSite");
							//System.out.println(readBaseOffsetHetSite);
							
							if(!(blockRefStart <= position && position <= (blockRefStart + blockLength -1))) continue;
							
							
							byte[] baseHetSite = record.getReadBases();
							//System.out.println(record.getReadString());
							//System.out.println((record.getReadString()).length());
							Character baseHet = new Character((char)baseHetSite[readBaseOffsetHetSite]);
							//System.out.println("here");
							
							for(String key : localCGSites.keySet()) {
								
								int pos = new Integer(key.split(",")[1]);
								//System.out.println(pos);
								if(bamStrand.equals("reverse")) {
									pos=pos+1;
								}
								//System.out.println(pos);
							
									//System.out.println(pos);
									//System.out.println(blockRefStart);
									//System.out.println(blockRefStart + blockLength -1);
								
								
								if(blockRefStart <= pos && pos <= (blockRefStart + blockLength -1)) { //block overlaps CG site
									//System.out.println(readName);
									CGSite cgSite = localCGSites.get(key);
									
									int refOffset = pos - blockRefStart;
									int readBaseOffset = blockReadStart -1 + refOffset;				
									Character baseChar = null;
									byte[] base = record.getReadBases();
									baseChar = new Character((char)base[readBaseOffset]);
									//System.out.println(baseChar);
									// Add if read is not already assigned to CG site
									//System.out.println(readName);
									//System.out.println(bamOrder);
									//System.out.println(bamStrand);
									if(!cgSite.checkRead(readName, bamOrder, bamStrand)) {
										//System.out.println(baseHet);
										if(baseHet == (char)cgSite.allele1.getBase()) {
											//System.out.println("here");
											addCount(cgSite, "allele1", bamOrder, bamStrand, baseChar);
									}
										else if(baseHet == (char)cgSite.allele2.getBase()) {
											//System.out.println("here");
											addCount(cgSite, "allele2", bamOrder, bamStrand, baseChar);
									}
										//System.out.println("here");
									cgSite.addRead(readName, bamOrder, bamStrand);
									}
									//System.out.println("here");
																		
									//System.out.println(pos);
									//System.out.println(baseHet);
									//System.out.println(baseChar);
									//byte baseHetSite = record.get
									//System.out.println(bamStrand);
									
								}
							}
						}
					}
					samIterator.close();
				}	
			}
			catch (RuntimeException e) {
				throw new RuntimeException("Problem reading file:" +e);
			} 
		}
		printCGSites(globalCGSites);
		
	}
	
	public Character invertAllele(Character allele) {
			if(allele.equals('A')) return 'T';
			else if(allele.equals('T')) return 'A';
			else if(allele.equals('G')) return 'C';
			else if(allele.equals('C')) return 'G';
			else return 'N';
	}
	
	public void addCount(CGSite cgsite, String allele, String bamOrder, String bamStrand, Character base) {
		
		if(allele.equals("allele1")) {
			if(bamOrder.equals("first")) {
				if(bamStrand.equals("forward")) {
					if(base == 'C') {
						cgsite.allele1.firstForwardCs++;
					} else if(base == 'T') {
						cgsite.allele1.firstForwardTs++;
					}
				} else if(bamStrand.equals("reverse")) {
					if(base == 'G') {
						cgsite.allele1.firstReverseCs++;
					} else if(base == 'A') {
						cgsite.allele1.firstReverseTs++;
					}
				}
			} else if(bamOrder.equals("second")) {
				if(bamStrand.equals("forward")) {
					if(base == 'C') {
						cgsite.allele1.secondForwardCs++;
					} else if(base == 'T') {
						cgsite.allele1.secondForwardTs++;
					}
				} else if(bamStrand.equals("reverse")) {
					if(base == 'G') {
						cgsite.allele1.secondReverseCs++;
					} else if(base == 'A') {
						cgsite.allele1.secondReverseTs++;
					}
				}
			}
		} else if(allele.equals("allele2")) {
			if(bamOrder.equals("first")) {
				if(bamStrand.equals("forward")) {
					if(base == 'C') {
						cgsite.allele2.firstForwardCs++;
					} else if(base == 'T') {
						cgsite.allele2.firstForwardTs++;
					}
				} else if(bamStrand.equals("reverse")) {
					if(base == 'G') {
						cgsite.allele2.firstReverseCs++;
					} else if(base == 'A') {
						cgsite.allele2.firstReverseTs++;
					}
				}
			} else if(bamOrder.equals("second")) {
				if(bamStrand.equals("forward")) {
					if(base == 'C') {
						cgsite.allele2.secondForwardCs++;
					} else if(base == 'T') {
						cgsite.allele2.secondForwardTs++;
					}
				} else if(bamStrand.equals("reverse")) {
					if(base == 'G') {
						cgsite.allele2.secondReverseCs++;
					} else if(base == 'A') {
						cgsite.allele2.secondReverseTs++;
					}
				}
			}
		}
		
	}
	
	public void printCGSites(HashMap<String, CGSite> globalCGSites) {
		
		System.out.println("chr,pos,a1,a2,cov_a1_first_forward,cov_a1_first_reverse,cov_a1_second_forward,cov_a1_second_reverse,cov_a2_first_forward," +
				"cov_a2_first_reverse,cov_a2_second_forward,cov_a2_second_reverse,met_a1_first_forward,met_a1_first_reverse,met_a1_second_forward," +
				"met_a1_second_reverse,met_a2_first_forward,met_a2_first_reverse,met_a2_second_forward,met_a2_second_reverse");
		
		
		for(String key : globalCGSites.keySet()) {
			String chr = key.split(",")[0];
			int pos = new Integer(key.split(",")[1]);
			CGSite cgsite = globalCGSites.get(key);
			Allele a1 = cgsite.allele1;
			Allele a2 = cgsite.allele2;

			double covFirstForwarda1 = a1.firstForwardTs+a1.firstForwardCs;
			double covFirstReversea1 = a1.firstReverseTs+a1.firstReverseCs;
			double covSecondForwarda1 = a1.secondForwardTs+a1.secondForwardCs;
			double covSecondReversea1 = a1.secondReverseTs+a1.secondReverseCs;
			double covFirstForwarda2 = a2.firstForwardTs+a2.firstForwardCs;
			double covFirstReversea2 = a2.firstReverseTs+a2.firstReverseCs;
			double covSecondForwarda2 = a2.secondForwardTs+a2.secondForwardCs;
			double covSecondReversea2 = a2.secondReverseTs+a2.secondReverseCs;
			
			double firstForwardMetha1 = -1;
			double firstReverseMetha1 = -1;
			double secondForwardMetha1 = -1;
			double secondReverseMetha1 = -1;
			double firstForwardMetha2 = -1;
			double firstReverseMetha2 = -1;
			double secondForwardMetha2 = -1;
			double secondReverseMetha2 = -1;
			
			if(covFirstForwarda1 != 0) firstForwardMetha1 = (a1.firstForwardCs/(covFirstForwarda1))*100;
			if(covFirstReversea1 != 0) firstReverseMetha1 = (a1.firstReverseCs/(covFirstReversea1))*100;
			if(covSecondForwarda1 != 0) secondReverseMetha1 = (a1.secondForwardCs/(covSecondForwarda1))*100;
			if(covSecondReversea1 != 0) secondReverseMetha1 = (a1.secondReverseCs/(covSecondReversea1))*100;
		
			if(covFirstForwarda2 != 0) firstForwardMetha2 = (a2.firstForwardCs/covFirstForwarda2)*100;
			if(covFirstReversea2 != 0) firstReverseMetha2 = (a2.firstReverseCs/covFirstReversea2)*100;
			if(covSecondForwarda2 != 0) secondForwardMetha2 = (a2.secondForwardCs/covSecondForwarda2)*100;
			if(covSecondReversea2 != 0) secondReverseMetha2 = (a2.secondReverseCs / covSecondReversea2)*100;

			firstForwardMetha1=roundNumber(firstForwardMetha1);
			firstReverseMetha1=roundNumber(firstReverseMetha1);
			secondReverseMetha1=roundNumber(secondReverseMetha1);
			secondReverseMetha1=roundNumber(secondReverseMetha1);
			firstForwardMetha2=roundNumber(firstForwardMetha2);
			firstReverseMetha2=roundNumber(firstReverseMetha2);
			secondForwardMetha2=roundNumber(secondForwardMetha2);
			secondReverseMetha2=roundNumber(secondReverseMetha2);
			
			 String toPrint = chr + "," + pos + "," + a1.getBase() + "," + a2.getBase() + "," + (int)covFirstForwarda1 + "," + (int)covFirstReversea1 + "," + 
			 + (int)covSecondForwarda1 + "," + (int)covSecondReversea1 + "," + (int)covFirstForwarda2 + "," + (int)covFirstReversea2 + "," + (int)covSecondForwarda2 + "," + (int)covSecondReversea2
			 + "," + firstForwardMetha1 + "," + firstReverseMetha1 + "," + secondForwardMetha1 + "," + secondReverseMetha1 + "," + firstForwardMetha2 + "," +
			firstReverseMetha2 + "," + secondForwardMetha2 + "," + secondReverseMetha2;
			
			/*String toPrint = chr + "," + pos + "," + a1.getBase() + "," + a2.getBase() + "," + covFirstForwarda1 + "," + covFirstReversea1 + "," 
			+ covSecondForwarda1 + "," + covSecondReversea1 + "," + covFirstForwarda2 + "," + covFirstReversea2 + "," + covSecondForwarda2 + "," + covSecondReversea2;
			*/
			System.out.println(toPrint);
		}
	}
	
	public double roundNumber(double a) {
		return Math.round(a * 100) / 100.0;
	}
}

			

			

