package ca.mcgill.genome.bamtools.commands;

import gnu.trove.map.TCharDoubleMap;
import gnu.trove.map.hash.TCharDoubleHashMap;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.parsers.CISBPParser;
import ca.mcgill.genome.bamtools.pwm.RnaBindingProteinInfo;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;

public class MotifCISBPRNA extends DefaultTool {
	private boolean printHeader = false;
	private double minQual = 0.0;
	private int sampleIndex = 0;
	private int minCLR = 0;
	private TCharDoubleMap atcgBGFrequencies = null;
	private IndexedFastaSequenceFile reference;
	private String sampleName = "";
	private SnpEffectPredictor predictor;

	@Override
	public String getCmdName() {
		return "motifCISBP";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("CISBP can be downloaded here: http://cisbp-rna.ccbr.utoronto.ca/");
		System.out.println("\t--vcf           VCF to test");
		System.out.println("\t--ref           Indexed Genome reference file");
		System.out.println("\t -c             SNPEff config file");
		System.out.println("\t--atcg          Comma separated genome background base frequencies. In A,T,C,G order. If not given, will be computed from --ref.");
		System.out.println("\t                Put 0.25,0.25,0.25,0.25 if you wish to skip this step.");
		System.out.println("\t--pwm           Position Weight Matrix directory for RBPs");
		System.out.println("\t--rbp           Rna Binding Protein info file");
		System.out.println("\t--bed           Rna Binding Protein genome coordinates (can be used multiple times)");
		System.out.println("\t--sampleIdx     Sample index in the VCF (default: " + sampleIndex + ")");
		System.out.println("\t--minQuality    Minimum variant quality (default: " + minQual + ")");
		System.out.println("\t--minCLR        Minimum CLR. (default: " + minCLR + ")");
		System.out.println("\t--header        Print the output header. Default doesn't print.");
		System.out.println("\t--sampleName    Sample Name to print");
	}

	@Override
	public int run(String[] args) {
		File vcfFile = null;
		List<File> beds = new ArrayList<File>();
		File pwmDir = null;
		File rbpFile = null;
		String atcgStr = null;
		File config = null;
		String snpEffGenome = null;

		try {
			for (int idx = 1; idx < args.length; idx++) {
				if (args[idx].equals("--vcf")) {
					idx++;
					vcfFile = new File(args[idx]);
				} else if (args[idx].equals("-c")) {
					idx++;
					config = new File(args[idx]);
					idx++;
					snpEffGenome = args[idx];
				} else if (args[idx].equals("--ref")) {
					idx++;
					reference = new IndexedFastaSequenceFile(new File(args[idx]));
				} else if (args[idx].equals("--atcg")) {
					idx++;
					atcgStr = args[idx];
				} else if (args[idx].equals("--pwm")) {
					idx++;
					pwmDir = new File(args[idx]);
				} else if (args[idx].equals("--rbp")) {
					idx++;
					rbpFile = new File(args[idx]);
				} else if (args[idx].equals("--bed")) {
					idx++;
					beds.add(new File(args[idx]));
				} else if (args[idx].equals("--minQuality")) {
					idx++;
					minQual = Double.parseDouble(args[idx]);
				} else if (args[idx].equals("--minCLR")) {
					idx++;
					minCLR = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--sampleIdx")) {
					idx++;
					sampleIndex = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--sampleName")) {
					idx++;
					sampleName = args[idx];
				} else if (args[idx].equals("--header")) {
					printHeader = true;
				}
				
			}

			if (vcfFile == null) {
				printUsage("vcf not set");
				return 1;
			}
			if (rbpFile == null) {
				printUsage("rbp not set");
				return 1;
			}
			if (beds.size() == 0) {
				printUsage("beds not set");
				return 1;
			}
			if (pwmDir == null) {
				printUsage("pwm not set");
				return 1;
			}
			if (!pwmDir.isDirectory()) {
				printUsage("pwm is not a directory");
				return 1;
			}

			Config snpEffConfig = new Config(snpEffGenome, config.toString());
			predictor =  snpEffConfig.loadSnpEffectPredictor();
			predictor.buildForest();

			if(atcgStr != null) {
				String atcgs[] = Gpr.split(atcgStr, ',');
				if(atcgs.length != 4) {
					printUsage("Need 4 values for atcg, A,T,C,G");
					return 1;
				}

				atcgBGFrequencies = new TCharDoubleHashMap();
				atcgBGFrequencies.put('A', Double.parseDouble(atcgs[0]));
				atcgBGFrequencies.put('T', Double.parseDouble(atcgs[1]));
				atcgBGFrequencies.put('C', Double.parseDouble(atcgs[2]));
				atcgBGFrequencies.put('G', Double.parseDouble(atcgs[3]));
			}
			else {
				atcgBGFrequencies = parseFrequencines();
			}
			System.err.print("Using frequencines A,T,C,G: ");
			System.err.print(atcgBGFrequencies.get('A'));
			System.err.print(',');
			System.err.print(atcgBGFrequencies.get('T'));
			System.err.print(',');
			System.err.print(atcgBGFrequencies.get('C'));
			System.err.print(',');
			System.err.print(atcgBGFrequencies.get('G'));
			System.err.println();
	
			CISBPParser parser = new CISBPParser(rbpFile, pwmDir);

			Map<String, RnaBindingProteinInfo> study2rbpsInformation  = parser.parse();
			IntervalForest forest = parseRegions(beds);
			testVariants(vcfFile, forest, study2rbpsInformation);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return 0;
	}

	private TCharDoubleMap parseFrequencines() {
		System.err.println("Computing BG frequencies");
		double totalnbBases = 0;
		double nbBases[] = new double[4];
		
		while(true) {
			ReferenceSequence refSeq = reference.nextSequence();
			if(refSeq == null)
				break;
			
			for(byte rawBaseValue : refSeq.getBases()) {
				char baseValue = Character.toUpperCase((char)rawBaseValue);
				
				// Count them independently to NOT count Ns.
				if(baseValue == 'A') {
					totalnbBases++;
					nbBases[0]++;
				}
				else if(baseValue == 'T') {
					totalnbBases++;
					nbBases[1]++;
				}
				else if(baseValue == 'C') {
					totalnbBases++;
					nbBases[2]++;
				}
				else if(baseValue == 'G') {
					totalnbBases++;
					nbBases[3]++;
				}
			}
		}
		TCharDoubleMap atcgFrequencies = new TCharDoubleHashMap();
		atcgFrequencies.put('A', nbBases[0]/totalnbBases);
		atcgFrequencies.put('T', nbBases[1]/totalnbBases);
		atcgFrequencies.put('C', nbBases[2]/totalnbBases);
		atcgFrequencies.put('G', nbBases[3]/totalnbBases);
		return atcgFrequencies;
	}

	private void printHeader() {
		System.out.print("SampleName");
		System.out.print("\t");
		System.out.print("Chromosome");
		System.out.print("\t");
		System.out.print("Start");
		System.out.print("\t");
		System.out.print("GeneName");
		System.out.print("\t");
		System.out.print("TargetGeneNames");
		System.out.print("\t");
		System.out.print("CISBPRNAId");
		System.out.print("\t");
		System.out.print("MotifId");
		System.out.print("\t");
		System.out.print("RawReferenceScore");
		System.out.print("\t");
		System.out.print("RawVariantScore");
		System.out.print("\t");
		System.out.print("RawLogRatio");
		System.out.print("\t");
		System.out.print("RelativeReferenceScore");
		System.out.print("\t");
		System.out.print("RelativeVariantScore");
		System.out.print("\t");
		System.out.print("RatioRelativeScores");
		System.out.println();
	}

	private void testVariants(File vcfFile, IntervalForest forest, Map<String, RnaBindingProteinInfo> study2rbpsInformation) {
		System.err.println("Testing variants");
		
		if(printHeader){
			printHeader();
		}
		
		VcfFileIterator vcfParser = new VcfFileIterator(vcfFile.toString());
		for (VcfEntry vcfEntry : vcfParser) {
			if(vcfEntry.getQuality() < minQual || vcfEntry.getInfoInt("CLR") < minCLR || !vcfEntry.isSnp()) {
				continue;
			}

			Markers intersects = forest.query(vcfEntry);
			if (intersects.size() > 0) {
				VcfGenotype sampleGenotype = vcfEntry.getVcfGenotype(sampleIndex);
				// test for hom ref
				if(sampleGenotype.getGenotypeCode() != 0) {
					for (Marker markerInt : intersects) {
						//????
						if(vcfEntry.getStart() < markerInt.getStart()) {
							continue;
						}
						RnaBindingProteinInfo rbpInfo = study2rbpsInformation.get(markerInt.getId());
						if(rbpInfo == null) {
							System.err.println("No RNA Binding Protein Info for: "+markerInt.getId());
							continue;
						}
//			System.err.println("Ref: " +vcfEntry.getChromosomeName()+':'+markerInt.getStart()+'-'+markerInt.getEnd()+" - " + vcfEntry.getStart());
						String refSequence;
						try {
							refSequence = new String(reference.getSubsequenceAt(vcfEntry.getChromosomeName(), markerInt.getStart()+1, markerInt.getEnd()+1).getBases(), "ASCII");
						} catch (UnsupportedEncodingException e) {
							throw new RuntimeException(e);
						}

//			System.err.println("Ref: "+refSequence);
						double refScore = rbpInfo.getPwm().score(refSequence);
						double refRelativeScore = rbpInfo.getPwm().relativeScore(refSequence, atcgBGFrequencies);
	
						int offset = vcfEntry.getStart() - markerInt.getStart();
						String variantSequence;
						if(offset == 0) {
//			System.err.println("Alt: 1");
							variantSequence = vcfEntry.getAlts()[0] + refSequence.substring(1);
						} else if(offset == rbpInfo.getPwm().size()) {
//			System.err.println("Alt: 2");
							variantSequence = refSequence.substring(0, offset) + vcfEntry.getAlts()[0];
						} else {
//			System.err.println("Alt: 3");
							variantSequence  = refSequence.substring(0, offset);
							variantSequence += vcfEntry.getAlts()[0];
							variantSequence += refSequence.substring(offset+1);
						}
//			System.err.println("Alt: "+variantSequence+'-'+vcfEntry.getRef()+'-'+vcfEntry.getAlts()[0]+'-'+offset+'-'+rbpInfo.getPwm().size());
						double varScore = rbpInfo.getPwm().score(variantSequence);
						double varRelativeScore = rbpInfo.getPwm().relativeScore(variantSequence, atcgBGFrequencies);
						
						Markers geneTargets = predictor.getIntervalForest().query(vcfEntry);						
						
						System.out.print(sampleName);
						System.out.print("\t");
						System.out.print(vcfEntry.getChromosomeName());
						System.out.print("\t");
						System.out.print(vcfEntry.getStart());
						System.out.print("\t");
						System.out.print(rbpInfo.getGeneName());
						System.out.print("\t");
						if(geneTargets.size() > 0) {
							Set<String> genes = new HashSet<String>();
							for(Marker marker : geneTargets) {
								if(marker instanceof ca.mcgill.mcb.pcingola.interval.Gene) {
									genes.add(((ca.mcgill.mcb.pcingola.interval.Gene)marker).getGeneName());
								}
							}
							boolean firstGene = true;
							for(String geneName : genes) {
								if(!firstGene) {
									System.out.print(',');
								}
								else {
									firstGene = false;
								}
								System.out.print(geneName);
							}
						}
						System.out.print("\t");
						
						
						System.out.print(rbpInfo.getCISBPRNAId());
						System.out.print("\t");
						System.out.print(rbpInfo.getMotifId());
						System.out.print("\t");
						System.out.print(refScore);
						System.out.print("\t");
						System.out.print(varScore);
						System.out.print("\t");
						System.out.print(Math.log(varScore/refScore));
						System.out.print("\t");
						System.out.print(refRelativeScore);
						System.out.print("\t");
						System.out.print(varRelativeScore);
						System.out.print("\t");
						System.out.print(varRelativeScore/refRelativeScore);
						System.out.println();
						
					}
				}
			}
		}
		
		vcfParser.close();
	}
	private IntervalForest parseRegions(List<File> regions) {
		System.err.println("Parsing Regions");
		IntervalForest forest = new IntervalForest();
		Genome genome = new Genome();

		for(File bed : regions) {
			LineNumberReader reader = null;
			try {
				reader = new LineNumberReader(new InputStreamReader(new FileInputStream(bed), "ASCII"));
				
				while(true) {
					String line = reader.readLine();
					if(line == null)
						break;
	
					// Ignore empty lines and comment lines
					if( (line.length() > 0) && (!line.startsWith("#")) ) {
						// Parse line
						String fields[] = line.split("\t");
	
						// Is line OK?
						if( fields.length >= 3 ) {
							// Format: CHR \t START \t END \t ID \t SCORE \t ....
							// Fields 
							String chromosome = fields[0].trim();
	
							// Start
							int start =  Gpr.parseIntSafe(fields[1]);
	
							// End
							int end = start; // Default 'end is same as start (only if value is missing)
							if( fields.length > 2 ) end = Gpr.parseIntSafe(fields[2]) - 1; // The chromEnd base is not included
	
							// ID
							String idFull = fields[3];
							int idx = idFull.indexOf('_');
							String id = idFull.substring(0, idx);
							if(id.length() < 8) {
								id = idFull.substring(0, idFull.indexOf('_', idx+1));
							}

							// Score and all following fields are ignored 

							// Create regulation
							forest.add(new Marker(genome.getOrCreateChromosome(chromosome), start, end, 1, id));
						}
					}
				}
			}
			catch(IOException e) {
				throw new RuntimeException(e);
			}
			finally {
				if(reader != null) {
					try {
						reader.close();
					} catch (IOException e) {
						//oh well.
					}
				}
			}
		}
		return forest;
	}
}
