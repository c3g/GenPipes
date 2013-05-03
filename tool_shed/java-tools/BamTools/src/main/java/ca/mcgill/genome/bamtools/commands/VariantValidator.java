package ca.mcgill.genome.bamtools.commands;

import gnu.trove.map.hash.TObjectDoubleHashMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.zip.GZIPInputStream;

import net.sf.picard.reference.IndexedFastaSequenceFile;

import org.apache.commons.lang.StringUtils;

import ca.mcgill.genome.bamtools.AlleleCounts;
import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.PairMappings;
import ca.mcgill.genome.bamtools.SampleAlleleCounts;
import ca.mcgill.genome.bamtools.SampleAlleleCounts.IndelAtPosition;
import ca.mcgill.genome.bamtools.Variant;
import ca.mcgill.genome.bamtools.parsers.ParseSimpleTSV;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;

public class VariantValidator extends DefaultTool {
	private double hetThreshold = 0.20;
	private int minDepth = 10;
	private double minQual = 0;
	private FormatVersion snpEffFormat = FormatVersion.FORMAT_SNPEFF_3;

	@Override
	public String getCmdName() {
		return "validator";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--vcfsFile           File containing a list of VCF files to get positions to validate from");
		System.out.println("\t--vcf                VCF to get positions to validate from");
		System.out.println("\t--regions            Target fragments to filter VCFs");
		System.out.println("\t--tsv                TSV file (sample, ref, alt, chr, pos)");
		System.out.println("\t--ref                Indexed reference genome");
		System.out.println("\t--validationbams     Samples BAM from the validation <needed with vtsv, tsv>");
		System.out.println("\t--referencebams      Samples BAM from another run <needed with vtsv, tsv>");
		System.out.println("\t--pairs              Pairs file <needed with vtsv>");
		System.out.println("\t--het                At which threshold do we call a valid het. (default: " + hetThreshold + ")");
		System.out.println("\t--mindepth           Minimum depth on which to call a Variant. (default: " + minDepth + ")");
		System.out.println("\t--minQual            Minimum variant quality to use to compute stats. (default: " + minQual + ")");
		System.out.println("\t--saveVariants       Minimum depth on which to call a Variant.");
		System.out.println("\t--loadVariants       Minimum depth on which to call a Variant.");
		System.out.println("\t--snpEffFormat       VCF snpEff format. (default: " + snpEffFormat + ")");
	}

	@Override
	public int run(String[] args) {
		try {
			List<File> vcfs = new ArrayList<File>();
			File vcfsFile = null;
			File regions = null;
			File tsv = null;
			IndexedFastaSequenceFile fastaRef = null;
			Map<String, File> validationSamplesBAM = null;
			Map<String, File> referenceSamplesBAM = null;
			PairMappings sampleMappings = null;
			File loadVariants = null; 
			File saveVariants = null; 
	
			for (int idx = 1; idx < args.length; idx++) {
				if (args[idx].equals("--tsv")) {
					idx++;
					tsv = new File(args[idx]);
				} else if (args[idx].equals("--vcf")) {
					idx++;
					vcfs.add(new File(args[idx]));
				} else if (args[idx].equals("--vcfsFile")) {
					idx++;
					vcfsFile = new File(args[idx]);
				} else if (args[idx].equals("--regions")) {
					idx++;
					regions = new File(args[idx]);
				} else if (args[idx].equals("--ref")) {
					idx++;
					fastaRef = new IndexedFastaSequenceFile(new File(args[idx]));
				} else if (args[idx].equals("--pairs")) {
					idx++;
					sampleMappings = PairMappings.parsePairsFile(new File(args[idx]));
				} else if (args[idx].equals("--validationbams")) {
					idx++;
					validationSamplesBAM = VariantAlleleFrequency.parseSamplesFile(new File(args[idx]));
				} else if (args[idx].equals("--referencebams")) {
					idx++;
					referenceSamplesBAM = VariantAlleleFrequency.parseSamplesFile(new File(args[idx]));
				} else if (args[idx].equals("--het")) {
					idx++;
					hetThreshold = Double.parseDouble(args[idx]);
				} else if (args[idx].equals("--snpEffFormat")) {
					idx++;
					snpEffFormat = Enum.valueOf(FormatVersion.class, args[idx]);
				} else if (args[idx].equals("--mindepth")) {
					idx++;
					minDepth = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--minQual")) {
					idx++;
					minQual = Double.parseDouble(args[idx]);
				} else if (args[idx].equals("--saveVariants")) {
					idx++;
					saveVariants = new File(args[idx]);
				} else if (args[idx].equals("--loadVariants")) {
					idx++;
					loadVariants = new File(args[idx]);
				}
			}
	
			if(vcfsFile != null) {
				BufferedReader reader = null;
				reader = new BufferedReader(new InputStreamReader(new FileInputStream(vcfsFile)));
				while(true) {
					String line = reader.readLine();
					if(line == null)
						break;
					if(line.startsWith("#"))
						continue;
					vcfs.add(new File(line));
				}
				reader.close();
			}
			
			if (vcfs.size() == 0 && tsv == null) {
				printUsage("vtsv and tsv not set");
				return 1;
			} else if (vcfs.size() > 0 && tsv != null) {
				printUsage("Both vtsv and tsv cannot be set");
				return 1;
			} else if (vcfs.size() > 0 && sampleMappings == null) {
				printUsage("pairs not set");
				return 1;
			} else if (vcfs.size() > 0 && regions == null) {
				printUsage("regions not set");
				return 1;
			} else if ((vcfs.size() > 0 || tsv != null) && validationSamplesBAM == null) {
				printUsage("bams not set");
				return 1;
			}
	
			List<Variant> variants=null;
			Map<String, SampleAlleleCounts> samplesBAMAlleleCounter = new HashMap<String, SampleAlleleCounts>();
			
			if(validationSamplesBAM != null) {
				for (String sampleName : validationSamplesBAM.keySet()) {
					File refSampleBAM = null;
					if(referenceSamplesBAM != null)
						refSampleBAM = referenceSamplesBAM.get(sampleName);
					samplesBAMAlleleCounter.put(sampleName, new SampleAlleleCounts(fastaRef, validationSamplesBAM.get(sampleName), refSampleBAM, 0, 0));
				}
			}
			
			if(loadVariants != null) {
				variants =  readVariants(loadVariants);
			}
			else if (vcfs.size() > 0) {
				variants = parseVariants(vcfs, regions, validationSamplesBAM.keySet());
			}
			else if (tsv != null) {
				variants = parseVariants(tsv);
			}
			
			if(saveVariants != null) {
				writeVariants(variants, saveVariants);
			}
			
			validateVariants(sampleMappings, hetThreshold, variants, samplesBAMAlleleCounter);
			for (SampleAlleleCounts counts : samplesBAMAlleleCounter.values()) {
				try {
					counts.close();
				} catch (IOException e) {
					// ignore
				}
			}

		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return 0;
	}

	private void writeVariants(List<Variant> variants, File saveVariants) {
		PrintWriter writer;
		try {
			writer = new PrintWriter(saveVariants, "ASCII");
			
			for(Variant variant : variants) {
				String dump = variant.write();
				writer.println(dump);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		
		writer.close();
	}

	private List<Variant> readVariants(File loadVariants) {
		List<Variant> variants = new ArrayList<Variant>();
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(loadVariants), "ASCII"));
			while(true) {
				String line = reader.readLine();
				if(line == null)
					break;
				if(line.startsWith("#"))
					continue;
				
				variants.add(Variant.read(line));
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		finally {
			try {
				reader.close();
			} catch (IOException e) {
				e.printStackTrace();
				//Do nothing
			}
		}
		return variants;
	}

	private IndelAtPosition getIndelAtPositionType(Variant variant) {
		IndelAtPosition indelAtPos = IndelAtPosition.NONE;
		int indelLength = variant.getRefAllele().length();
		// If it is an indel give the ref allele length
		for(String altAllele : variant.getAltAlleles()) {
			if(altAllele.length() > indelLength) {
				if(indelAtPos == IndelAtPosition.NONE) {
					indelAtPos = IndelAtPosition.INS;
				}
				else if(indelAtPos == IndelAtPosition.DEL) {
					indelAtPos = IndelAtPosition.BOTH;
				}
			}
			else if(altAllele.length() < indelLength) {
				if(indelAtPos == IndelAtPosition.NONE) {
					indelAtPos = IndelAtPosition.DEL;
				}
				else if(indelAtPos == IndelAtPosition.INS) {
					indelAtPos = IndelAtPosition.BOTH;
				}
			}
			else if(indelLength > 1 && altAllele.length() == indelLength) {
				throw new RuntimeException("MNPs are not handled");
			}
		}
		return indelAtPos;
	}

	private void validateVariants(PairMappings sampleMapping, double hetThreshold, List<Variant> variants, Map<String, SampleAlleleCounts> samplesBAMAlleleCounter) {
		long totalTests = 0;
		long totalSampleVariantTests = 0;
		long totalPassed = 0;
		long totalSampleVariantPassed = 0;
		long totalNonIndelTests = 0;
		long totalNonIndelPassed = 0;
		long totalFP = 0;
		long totalFN = 0;
		System.out.println("Sample\tChr\tPos\tGene\tVariant\tSampleAlleles\tOriginalQual\tOriginalCLR\tReferenceBAMAlleleCounts\tAlleleCounts\tReferenceAlleleFractions\tAlleleFractions\tVerdict\tMessage");
		for (Variant variant : variants) {
			for (Variant.SampleGenotype sampleGenotype : variant.getSampleGenotypes()) {
				SampleAlleleCounts sampleAlleleCounts = samplesBAMAlleleCounter.get(sampleGenotype.getSampleName());
				if(sampleAlleleCounts == null)
					continue;
				
				IndelAtPosition indelAtPos = getIndelAtPositionType(variant);
				int indelLength = variant.getRefAllele().length();

				Callable<AlleleCounts> callableCount = sampleAlleleCounts.getCallableAlleleCountsComputer(variant.getChromosome(), variant.getPosition(), indelAtPos, indelLength);
				AlleleCounts baseCounts;
				try {
					baseCounts = callableCount.call();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
				StringBuilder valSB = baseCounts.format(false);
				StringBuilder refSB;

				TObjectDoubleHashMap<String> valFraction = AlleleCounts.computeAlleleFractions(baseCounts.getBamCounts());
				StringBuilder valFracSB = AlleleCounts.formatFractionCounts(valFraction);
				StringBuilder refFracSB;

				if(baseCounts.getReferenceBAMCounts() != null) {
					refSB = baseCounts.format(true);
					TObjectDoubleHashMap<String> refFraction = AlleleCounts.computeAlleleFractions(baseCounts.getReferenceBAMCounts());
					refFracSB = AlleleCounts.formatFractionCounts(refFraction);
				}
				else {
					refSB = new StringBuilder();
					refFracSB = new StringBuilder();
				}

				int depth=baseCounts.getDepth();

				if(sampleMapping != null) {
					String sampleName = sampleMapping.getSampleFromNormal(sampleGenotype.getSampleName());
					if(sampleName != null) {
						System.out.print(sampleName+"-N");
					}
					else {
						sampleName = sampleMapping.getSampleFromTumor(sampleGenotype.getSampleName());
						if(sampleName != null) {
							System.out.print(sampleName+"-T");
						}
						else {
							System.out.print(sampleGenotype.getSampleName());
						}
					}
				}
				else {
					System.out.print(sampleGenotype.getSampleName());
				}
				System.out.print('\t');
				System.out.print(variant.getChromosome());
				System.out.print('\t');
				System.out.print(variant.getPosition());
				System.out.print('\t');
				System.out.print(variant.getGene());
				System.out.print('\t');
				System.out.print(variant.getRefAllele());
				System.out.print(" / ");
				System.out.print(StringUtils.join(variant.getAltAlleles(), ','));
				System.out.print('\t');
				
				if(sampleGenotype.getAlleleA() == null) {
					System.out.print("NA");
				}
				else {
					System.out.print(sampleGenotype.getAlleleA());
					System.out.print(" / ");
					System.out.print(sampleGenotype.getAlleleB());
				}
				System.out.print('\t');
				System.out.print(sampleGenotype.getQual());
				System.out.print('\t');
				System.out.print(sampleGenotype.getClr());
				System.out.print('\t');
				System.out.print(refSB);
				System.out.print('\t');
				System.out.print(valSB);
				System.out.print('\t');
				System.out.print(refFracSB);
				System.out.print('\t');
				System.out.print(valFracSB);
				System.out.print('\t');

				if(sampleGenotype.getAlleleA() == null) {
					TestResult tstResult = testVariant(hetThreshold, variant, variant.getRefAllele(), variant.getAltAlleles().get(0), valFraction);
					if(tstResult.getOutcome() == TestResult.TestOutcome.VALID) {
						System.out.print("Would have a HET variant");						
					}
					else {
						tstResult = testVariant(hetThreshold, variant, variant.getAltAlleles().get(0), variant.getAltAlleles().get(0), valFraction);
						if(tstResult.getOutcome() == TestResult.TestOutcome.VALID) {
							System.out.print("Would have a HOM variant");						
						}
						else {
							System.out.print("Might not have a variant");
						}
					}
					System.out.print('\t');
					System.out.print("Not previously sequenced");
				}
				else if(depth >= minDepth) {
					totalTests++;
					if (!variant.isIndel()) {
						totalNonIndelTests++;
					}

					if(sampleGenotype.isVariant(variant.getRefAllele()) && sampleGenotype.getQual() >= minQual) {
						totalSampleVariantTests++;
					}

					TestResult tstResult = testVariant(hetThreshold, variant, sampleGenotype.getAlleleA(), sampleGenotype.getAlleleB(), valFraction);
					if (tstResult.getOutcome() == TestResult.TestOutcome.VALID) {
						System.out.print("PASSED");
						System.out.print('\t');
						totalPassed++;
						if (!variant.isIndel()) {
							totalNonIndelPassed++;
						}
						if(sampleGenotype.isVariant(variant.getRefAllele()) && sampleGenotype.getQual() >= minQual) {
							totalSampleVariantPassed++;
						}
					} else if (tstResult.getOutcome() == TestResult.TestOutcome.FP) {
						System.out.print("FAILED-FP");
						System.out.print('\t');
						System.out.print(tstResult.getMessage());
						if (!variant.isIndel()) {
							totalFP++;
						}
					} else if (tstResult.getOutcome() == TestResult.TestOutcome.FN) {
						System.out.print("FAILED-FN");
						System.out.print('\t');
						System.out.print(tstResult.getMessage());
						if (!variant.isIndel()) {
							totalFN++;
						}
					} else if (tstResult.getOutcome() == TestResult.TestOutcome.NA) {
						System.out.print("Unknown outcome");
						System.out.print('\t');
						System.out.print(tstResult.getMessage());
					} else {
						throw new RuntimeException("Unknown outcome: " + tstResult.getOutcome());
					}
				}
				else {
					System.out.print("IGNORED");
					System.out.print('\t');
					System.out.print("Depth too low");
				}
				System.out.println();
			}
		}
		DecimalFormat formatter = new DecimalFormat("#0.00");
		System.out.println("Sample Variant Passed:  " + formatter.format((double)totalSampleVariantPassed * 100d / (double)totalSampleVariantTests) + "%");
		System.out.println("Passed:                 " + formatter.format((double)totalPassed * 100d / (double)totalTests) + "%");
		System.out.println("Non-indel Passed:       " + formatter.format((double)totalNonIndelPassed * 100d / (double)totalNonIndelTests) + "%");
		System.out.println("FP percent:             " + formatter.format((double)totalFP * 100d / (double)totalNonIndelTests) + "%");
		System.out.println("FN percent:             " + formatter.format((double)totalFN * 100d / (double)totalNonIndelTests) + "%");
	}

	private TestResult testVariant(double hetThreshold, Variant variant, String alleleA, String alleleB, TObjectDoubleHashMap<String> baseFraction) {
		TestResult.TestOutcome badOutcome;
		if(variant.getRefAllele().equals(alleleA) && variant.getRefAllele().equals(alleleB)) {
			badOutcome = TestResult.TestOutcome.FN;
		}
		else {
			badOutcome = TestResult.TestOutcome.FP;
		}

		TestResult result = new TestResult();
		result.setOutcome(TestResult.TestOutcome.VALID);
		List<String> keptBases = new ArrayList<String>();
		for (String base : baseFraction.keySet()) {
			if (baseFraction.get(base) > hetThreshold) {
				if (keptBases.size() == 2) {
					result.setMessage("Too many bases above het threshold [" + hetThreshold + "]");
					result.setOutcome(TestResult.TestOutcome.NA);
				} else {
					keptBases.add(base.toString());
				}
			}
		}

		if (!keptBases.contains(alleleA)) {
			result.setMessage("Missing allele '" + alleleA + "'");
			result.setOutcome(badOutcome);
		}
		if (!keptBases.contains(alleleB)) {
			result.setMessage("Missing allele '" + alleleB + "'");
			result.setOutcome(badOutcome);
		}

		keptBases.remove(alleleA);
		keptBases.remove(alleleB);

		if (keptBases.size() != 0) {
			result.setMessage("Spurrious allele " + StringUtils.join(keptBases, '/'));
			result.setOutcome(badOutcome);
		}
		return result;
	}

	private IntervalForest parseRegions(File regions) {
		IntervalForest forest = new IntervalForest();
		LineNumberReader reader = null;
		Genome genome = new Genome();
		try {
			reader = new LineNumberReader(new InputStreamReader(new FileInputStream(regions), "ASCII"));
			
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
						String id = "line_" + reader.getLineNumber();
						if( (fields.length > 3) && (!fields[3].isEmpty()) ) id = fields[3];

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
		return forest;
	}

	private List<Variant> parseVariants(List<File> inputVCFs, File regions, Set<String> samplesWithBAM) {
		Map<String, Variant> retVal = new HashMap<String, Variant>();

		IntervalForest forestOfRegions = parseRegions(regions);
		
		Set<String> samplesWithGT = new HashSet<String>();
		try {
			for(File vcf : inputVCFs) {
				System.err.println("Reading: " + vcf);
				String vcfName = vcf.toString();
				BufferedReader in;
				if(vcfName.endsWith(".gz")) {
					in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(vcf), 1024*1024), "ASCII"));
				}
				else {
					in = new BufferedReader(new InputStreamReader(new FileInputStream(vcf), "ASCII"), 1024*1024);
				}
				VcfFileIterator vcfFile = new VcfFileIterator(in);
				VcfHeader header = vcfFile.readHeader();
				samplesWithGT.add(header.getSampleNames().get(0));
				samplesWithGT.add(header.getSampleNames().get(1));
				for(VcfEntry vcfEntry : vcfFile) {
					String key = vcfEntry.getChromosomeName()+':'+(vcfEntry.getStart()+1);
					if(retVal.containsKey(key)) {
						Variant variant = retVal.get(key);
						if(!variant.getAltAlleles().get(0).equals(vcfEntry.getAlts()[0])) {
							System.err.println("Different alt alleles: "+key + " alt1: " + variant.getAltAlleles().get(0) + ' ' + vcfEntry.getAltsStr());
						}
						VcfGenotype normGT = vcfEntry.getVcfGenotypes().get(0);
						VcfGenotype tumorGT = vcfEntry.getVcfGenotypes().get(1);
						int clr = (int)vcfEntry.getInfoInt("CLR");
						double qual = vcfEntry.getQuality();

						variant.addSampleGenotype(header.getSampleNames().get(0), normGT.getGenotype(0), normGT.getGenotype(1), qual, clr);
						variant.addSampleGenotype(header.getSampleNames().get(1), tumorGT.getGenotype(0), tumorGT.getGenotype(1), qual, clr);
					}
					else {
						Markers hits = forestOfRegions.query(vcfEntry);
						if(hits.size() > 0) {
							List<VcfEffect> effects = vcfEntry.parseEffects(snpEffFormat);
							String gene = null;
							for(VcfEffect effect : effects) {
								if(effect.getGene() != null) {
									gene = effect.getGene();
									break;
								}
							}
							int clr = (int)vcfEntry.getInfoInt("CLR");
							double qual = vcfEntry.getQuality();

							Variant variant = new Variant(vcfEntry.getChromosomeName(), vcfEntry.getStart()+1, vcfEntry.getRef(), vcfEntry.getAltsStr(), gene);
							variant.simplifyIndels();
							retVal.put(key, variant);
							
							VcfGenotype normGT = vcfEntry.getVcfGenotypes().get(0);
							VcfGenotype tumorGT = vcfEntry.getVcfGenotypes().get(1);
							variant.addSampleGenotype(header.getSampleNames().get(0), normGT.getGenotype(0), normGT.getGenotype(1), qual, clr);
							variant.addSampleGenotype(header.getSampleNames().get(1), tumorGT.getGenotype(0), tumorGT.getGenotype(1), qual, clr);
						}
					}
				}
				vcfFile.close();
			}
		}
		catch(IOException e) {
			throw new RuntimeException(e);
		}
		
		List<String> notPrevCalled = new ArrayList<String>();
		for(String sampleWithBAM : samplesWithBAM){
			if(!samplesWithGT.contains(sampleWithBAM)) {
				notPrevCalled.add(sampleWithBAM);
			}
		}

		List<Variant> variants = new ArrayList<Variant>(retVal.values());
		for(Variant variant : variants) {
			for(String sample : samplesWithGT) {
				if(variant.getSampleToGenotype().get(sample) == null) {
					variant.addSampleGenotype(sample, variant.getRefAllele(), variant.getRefAllele(), -1, -1);
				}
			}
			for(String sample : notPrevCalled) {
				variant.addSampleGenotype(sample, null, null, -1, -1);
			}
		}
		
		return variants;
	}

	private List<Variant> parseVariants(File inputTSV) {
		List<Variant> retVal = new ArrayList<Variant>();

		ParseSimpleTSV parseVariantTSV = null;
		try {
			parseVariantTSV = new ParseSimpleTSV(inputTSV);
			parseVariantTSV.parseHeader();
			String nextLine[];
			try {
				while ((nextLine = parseVariantTSV.readNext()) != null) {
					Variant variant = new Variant(nextLine[parseVariantTSV.getChromosomeIdx()], Integer.parseInt(nextLine[parseVariantTSV.getPositionIdx()].replace(",", "")), nextLine[parseVariantTSV.getRefAlleleIdx()], nextLine[parseVariantTSV.getAltAllelesIdx()]);
					variant.simplifyIndels();
					retVal.add(variant);
					variant.addSampleGenotype(nextLine[parseVariantTSV.getSampleIdx()], variant.getRefAllele(), variant.getAltAlleles().get(0), -1, -1);
				}
			} catch (NumberFormatException e) {
				throw new RuntimeException(e);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		} finally {
			try {
				parseVariantTSV.close();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}
		return retVal;
	}
	
	private static class TestResult {
		public enum TestOutcome {FP, FN, VALID, NA};

		private TestOutcome outcome;
		private String message;
		public TestOutcome getOutcome() {
			return outcome;
		}
		public void setOutcome(TestOutcome outcome) {
			this.outcome = outcome;
		}
		public String getMessage() {
			return message;
		}
		public void setMessage(String message) {
			this.message = message;
		}
	}
}
