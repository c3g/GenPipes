package ca.mcgill.genome.bamtools.commands;

import gnu.trove.map.hash.TObjectLongHashMap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.mutations.MutationRateTarget;
import ca.mcgill.genome.bamtools.mutations.Region;
import ca.mcgill.genome.bamtools.parsers.GeneMapper;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Cds;
import ca.mcgill.mcb.pcingola.interval.Chromosome;
import ca.mcgill.mcb.pcingola.interval.Exon;
import ca.mcgill.mcb.pcingola.interval.Gene;
import ca.mcgill.mcb.pcingola.interval.Genes;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Intergenic;
import ca.mcgill.mcb.pcingola.interval.Intron;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.Transcript;
import ca.mcgill.mcb.pcingola.interval.Utr;
import ca.mcgill.mcb.pcingola.interval.Utr3prime;
import ca.mcgill.mcb.pcingola.interval.Utr5prime;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalTree;
import ca.mcgill.mcb.pcingola.snpEffect.ChangeEffect.EffectType;
import ca.mcgill.mcb.pcingola.snpEffect.Config;
import ca.mcgill.mcb.pcingola.snpEffect.SnpEffectPredictor;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;

public class MutationRates extends DefaultTool {
	public static final int NORMAL_VCF_IDX = 0;
	public static final int TUMOR_VCF_IDX = 1;
	private List<File> vcfFiles;
	private File output;
	private SnpEffectPredictor predictor;
	private int minGenotypeQual = 0;
	private double minQual = 0.0;
	private int minCLR = 45;
//	private final List<Region> regions = new ArrayList<Region>();
	private Region dnaseRegion;
	private long totalSize;

	@Override
	public String getCmdName() {
		return "mutationrates";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--vcfsFile     File containing a list of VCF files to get positions to validate from");
		System.out.println("\t--vcf          VCFs to get regions for");
		System.out.println("\t -c            SNPEff config file");
		System.out.println("\t--dnase        DNase regions");
		System.out.println("\t--mapper       Genelist to use (filter by refseq)");
		System.out.println("\t--minQual      Minimum variant quality. (Default: "+minQual+")");
		System.out.println("\t--minGQ        Minimum genotype quality. (Default: "+minGenotypeQual+")");
		System.out.println("\t--minCLR       Minimum CLR to call somatic. (Default: "+minCLR+")");
		System.out.println("\t--totalSize    Genome size (optional)");
		System.out.println("\t--output       Output file for stats");
		//System.out.println("\t--printROI   Dump RegionsOfInterest");
	}

	@Override
	public int run(String[] args) {
		File vcfsFile = null;
		vcfFiles = new ArrayList<File>();
		File config = null;
		String snpEffGenome = null;
		GeneMapper mapper = null;
		totalSize=(long)(3000000000d*0.9);

		try {
			for (int idx = 0; idx < args.length; idx++) {
				if (args[idx].equals("--vcf")) {
					idx++;
					vcfFiles.add(new File(args[idx]));
				} else if (args[idx].equals("--vcfsFile")) {
					idx++;
					vcfsFile = new File(args[idx]);
				} else if (args[idx].equals("--output")) {
					idx++;
					output = new File(args[idx]);
				} else if (args[idx].equals("-c")) {
					idx++;
					config = new File(args[idx]);
					idx++;
					snpEffGenome = args[idx];
				} else if (args[idx].equals("--dnase")) {
					idx++;
					dnaseRegion = parseDNAse(new File(args[idx]));
				} else if (args[idx].equals("--mapper")) {
					idx++;
					mapper = new GeneMapper(new File(args[idx]));
				} else if (args[idx].equals("--minGQ")) {
					idx++;
					minGenotypeQual = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--minQual")) {
					idx++;
					minQual = Double.parseDouble(args[idx]);
				} else if (args[idx].equals("--minCLR")) {
					idx++;
					minCLR = Integer.parseInt(args[idx]);
				} else if (args[idx].equals("--totalSize")) {
					idx++;
					totalSize = Long.parseLong(args[idx]);
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
					vcfFiles.add(new File(line));
				}
				reader.close();
			}
	
			if(vcfFiles.size() == 0) {
				printUsage("Missing vcfFiles");
				return 1;
			}
			if(output == null) {
				printUsage("Missing output");
				return 1;
			}
	
			Config snpEffConfig = new Config(snpEffGenome, config.toString());
			predictor =  snpEffConfig.loadSnpEffectPredictor();
			Genome genome = predictor.getGenome();
			
			Genes genes = genome.getGenes();
			if(mapper != null) {
				for(Gene gene : genes) {
					Iterator<Transcript> trIter = gene.iterator();
					while(trIter.hasNext()) {
						Transcript tr = trIter.next();
						if(mapper.getTranscriptId(tr.getId()) == null) {
							trIter.remove();
						}
						else if(!tr.getBioType().equalsIgnoreCase("protein_coding") && !tr.getBioType().equalsIgnoreCase("mRNA")) {
							trIter.remove();
						}
					}
				}
			}
	
			Iterator<Gene> geneIter = genes.iterator();
			Map<String, Integer> ends = new HashMap<String, Integer>();
			while(geneIter.hasNext()) {
				Gene gene = geneIter.next();
				if(gene.numChilds() == 0) {
					geneIter.remove();
				}
				else {
					List<Transcript> trs = gene.sorted();
					if(trs.get(0).getStart() > gene.getStart()) {
						gene.setStart(trs.get(0).getStart());
					}
					if(trs.get(trs.size()-1).getEnd() < gene.getEnd()) {
						gene.setEnd(trs.get(trs.size()-1).getEnd());
					}
	
					if(!ends.containsKey(gene.getChromosomeName())) {
						ends.put(gene.getChromosomeName(), gene.getEnd());
					}
					if(ends.get(gene.getChromosomeName()).intValue() < gene.getEnd()) {
						ends.put(gene.getChromosomeName(),gene.getEnd());
					}
				}
			}
	
			for(String chrName : ends.keySet()) {
				int endPos = ends.get(chrName);
				predictor.add(new Intergenic(genome.getChromosome(chrName), endPos, genome.getChromosome(chrName).getEnd(), 1, "..."));
			}
			predictor.buildForest();
	
			computeRates();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return 0;
	}
	
	public void computeRegionSizes(TObjectLongHashMap<String>  regionSizes) {
		regionSizes.put(dnaseRegion.getName(), dnaseRegion.getRegionSize());
		
		Markers utr5Prime = new Markers();
		Markers utr3Prime = new Markers();
		Markers cds = new Markers();
		Markers transcripts = new Markers();
		Markers upstream = new Markers();
		Markers downstream = new Markers();
		Markers introns = new Markers();
		Markers exons = new Markers();

		for(Gene gene : predictor.getGenome().getGenes()) {
			for (Transcript tr : gene) {
				int intronLength = 0;
				int exonLength = 0;
				
				for (Utr5prime utr : tr.get5primeUtrs())
					utr5Prime.add(utr);
				for (Utr3prime utr : tr.get3primeUtrs())
					utr3Prime.add(utr);
				for (Cds trCds : tr.getCds())
					cds.add(trCds);
				for (Intron intron : tr.introns()) {
					introns.add(intron);
					intronLength += intron.size();
				}
				for (Exon exon : tr) {
					exons.add(exon);
					exonLength += exon.size();
				}
				int estIntron = tr.size() - exonLength;
				if(estIntron != intronLength) {
					System.err.println("Not a good intron estimator");
				}
				transcripts.add(tr);
				upstream.add(tr.getUpstream());
				downstream.add(tr.getDownstream());
			}
			if(gene.sizeof(EffectType.INTRAGENIC.toString()) > 0) {
				System.err.println("Gene has intragenic: " + gene.getGeneName());
			}
		}
		
		Markers utr5PrimeMerged = utr5Prime.merge();
		Markers utr3PrimeMerged = utr3Prime.merge();
		Markers cdsMerged = cds.merge();
		Markers intronsMerged = introns.merge();
		//Markers transcriptsMerged = transcripts.merge();
		Markers upstreamMerged = upstream.merge();
		Markers downstreamMerged = downstream.merge();

		int size;
		size=0;
		for (Marker marker : utr5PrimeMerged) {
			size += marker.size();
		}
		regionSizes.adjustOrPutValue(EffectType.UTR_5_PRIME.toString(), size, size);
		size=0;
		for (Marker marker : utr3PrimeMerged) {
			size += marker.size();
		}
		regionSizes.adjustOrPutValue(EffectType.UTR_3_PRIME.toString(), size, size);
		size=0;
		for (Marker marker : cdsMerged) {
			size += marker.size();
		}
		regionSizes.adjustOrPutValue(EffectType.CDS.toString(), size, size);
		size=0;
		for (Marker marker : intronsMerged) {
			size += marker.size();
		}
		regionSizes.adjustOrPutValue(EffectType.INTRON.toString(), size, size);
		size=0;
		for (Marker marker : upstreamMerged) {
			size += marker.size();
		}
		regionSizes.adjustOrPutValue(EffectType.UPSTREAM.toString(), size, size);
		size=0;
		for (Marker marker : downstreamMerged) {
			size += marker.size();
		}
		regionSizes.adjustOrPutValue(EffectType.DOWNSTREAM.toString(), size, size);
	}

	public void computeRates() {
		MutationRateTarget somaticMutRateTarget = new MutationRateTarget();
		MutationRateTarget tumorMutRateTarget = new MutationRateTarget();
		MutationRateTarget normalMutRateTarget = new MutationRateTarget();

		TObjectLongHashMap<String> regionSizes = new TObjectLongHashMap<String>();
		computeRegionSizes(regionSizes);

		String type = EffectType.INTERGENIC.toString();
		for(IntervalTree tree : predictor.getIntervalForest()) {
			for(Marker marker : tree) {
				if(!marker.getChromosomeName().startsWith("H") && !marker.getChromosomeName().startsWith("G") && marker.getType() == EffectType.INTERGENIC || marker.getType() == EffectType.INTERGENIC_CONSERVED) {
					regionSizes.adjustOrPutValue(type, marker.size(), marker.size());
				}
			}
		}

		for(File vcfFile : vcfFiles) {
			String sampleName = vcfFile.getName().substring(0, vcfFile.getName().indexOf('.'));
			System.out.println("Processing " + sampleName + "...");
			
			somaticMutRateTarget.getCountPerTypePerSample().put(sampleName, new TObjectLongHashMap<String>());
			somaticMutRateTarget.getCountPerTypePerSample().get(sampleName).put("total", 0l);
			tumorMutRateTarget.getCountPerTypePerSample().put(sampleName, new TObjectLongHashMap<String>());
			tumorMutRateTarget.getCountPerTypePerSample().get(sampleName).put("total", 0l);
			normalMutRateTarget.getCountPerTypePerSample().put(sampleName, new TObjectLongHashMap<String>());
			normalMutRateTarget.getCountPerTypePerSample().get(sampleName).put("total", 0l);
			
			VcfFileIterator vcfParser = new VcfFileIterator(vcfFile.toString());
			TObjectLongHashMap<String> somaticSubstitutionCounts = new TObjectLongHashMap<String>();
			TObjectLongHashMap<String> tumorSubstitutionCounts = new TObjectLongHashMap<String>();
			TObjectLongHashMap<String> normalSubstitutionCounts = new TObjectLongHashMap<String>();
			for(String var : MutationRateTarget.substitutionOrder) {
				somaticSubstitutionCounts.put(var, 0);
				tumorSubstitutionCounts.put(var, 0);
				normalSubstitutionCounts.put(var, 0);
			}
			somaticMutRateTarget.getSubstitutionCountPerSample().put(sampleName, somaticSubstitutionCounts);
			tumorMutRateTarget.getSubstitutionCountPerSample().put(sampleName, tumorSubstitutionCounts);
			normalMutRateTarget.getSubstitutionCountPerSample().put(sampleName, normalSubstitutionCounts);
			somaticMutRateTarget.getSubstitutionCountPerTypePerSample().put(sampleName, new TObjectLongHashMap<String>());
			tumorMutRateTarget.getSubstitutionCountPerTypePerSample().put(sampleName, new TObjectLongHashMap<String>());
			normalMutRateTarget.getSubstitutionCountPerTypePerSample().put(sampleName, new TObjectLongHashMap<String>());
			
			try {
				for (VcfEntry vcfEntry : vcfParser) {
					if(vcfEntry.getQuality() < minQual) {
						continue;
					}

					VcfGenotype normalGenotype = vcfEntry.getVcfGenotype(NORMAL_VCF_IDX);
					VcfGenotype tumorGenotype = vcfEntry.getVcfGenotype(TUMOR_VCF_IDX);
					boolean isSomatic = vcfEntry.getInfoInt("CLR") >= minCLR;
					if(vcfEntry.isSnp()) {
						String substKey = vcfEntry.getRef()+'>'+vcfEntry.getAlts()[0];

						boolean normalHasVariant = false;
						if(Gpr.parseIntSafe(normalGenotype.get("GQ")) >= minGenotypeQual && normalGenotype.getGenotypeCode() != 0) {
							normalHasVariant = true;
							if(normalMutRateTarget.getSubstitutionCountPerSample().get(sampleName).increment(substKey) == false) {
								throw new RuntimeException("Missing titV" + vcfEntry);
							}
							normalMutRateTarget.getSubstitutionCountPerSample().get(sampleName).increment("total");
						}

						boolean tumorHasVariant = false;
						if(Gpr.parseIntSafe(tumorGenotype.get("GQ")) >= minGenotypeQual && tumorGenotype.getGenotypeCode() != 0) {
							tumorHasVariant = true;
							if(tumorMutRateTarget.getSubstitutionCountPerSample().get(sampleName).increment(substKey) == false) {
								throw new RuntimeException("Missing titV" + vcfEntry);
							}
							tumorMutRateTarget.getSubstitutionCountPerSample().get(sampleName).increment("total");
						}
						else {
							isSomatic = false;
						}
						
						if(isSomatic) {
							if(somaticMutRateTarget.getSubstitutionCountPerSample().get(sampleName).increment(substKey) == false) {
								throw new RuntimeException("Missing titV" + vcfEntry);
							}
							somaticMutRateTarget.getSubstitutionCountPerSample().get(sampleName).increment("total");

						}

						if(dnaseRegion.intersects(vcfEntry)) {
							String dnaseSubstKey = substKey+'_'+dnaseRegion.getName();
							
							if(isSomatic) {
								somaticMutRateTarget.getCountPerTypePerSample().get(sampleName).adjustOrPutValue(dnaseRegion.getName(), 1, 1);
								somaticMutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).adjustOrPutValue(dnaseSubstKey, 1, 1);
							}
							if(tumorHasVariant) {
								tumorMutRateTarget.getCountPerTypePerSample().get(sampleName).adjustOrPutValue(dnaseRegion.getName(), 1, 1);
								tumorMutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).adjustOrPutValue(dnaseSubstKey, 1, 1);
							}
							if(normalHasVariant) {
								normalMutRateTarget.getCountPerTypePerSample().get(sampleName).adjustOrPutValue(dnaseRegion.getName(), 1, 1);
								normalMutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).adjustOrPutValue(dnaseSubstKey, 1, 1);
							}
						}

						Set<Marker> markers = getHits(vcfEntry);
						if(markers.size() == 0) {
							System.err.println("Strange: "+vcfEntry);
						}

						Set<String> regionCounted = new HashSet<String>();
						for(Marker marker : markers) {
							String typeName = marker.getType().toString();
							// We only ignore splice site types since we test utr, cds and introns, we should have all the others.
							if(!regionSizes.containsKey(typeName) && !typeName.startsWith("SPLICE_SITE_")) {
								System.err.println("Doesn't contain type: " + typeName);
							}

							// Don't count the same region type more than once.
							if(regionCounted.contains(typeName))
								continue;
							regionCounted.add(typeName);
							String regionSubstKey = substKey+'_'+typeName;
							
							if(isSomatic) {
								somaticMutRateTarget.getCountPerTypePerSample().get(sampleName).adjustOrPutValue(typeName, 1, 1);
								somaticMutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).adjustOrPutValue(regionSubstKey, 1, 1);
							}
							if(tumorHasVariant) {
								tumorMutRateTarget.getCountPerTypePerSample().get(sampleName).adjustOrPutValue(typeName, 1, 1);
								tumorMutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).adjustOrPutValue(regionSubstKey, 1, 1);
							}
							if(normalHasVariant) {
								normalMutRateTarget.getCountPerTypePerSample().get(sampleName).adjustOrPutValue(typeName, 1, 1);
								normalMutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).adjustOrPutValue(regionSubstKey, 1, 1);
							}
						}
					}
				}
			}	
			finally {
				vcfParser.close();
			}
		}//for samples
		
		List<String> orderedNames = new ArrayList<String>();
		for(String regionName : regionSizes.keySet()) {
			orderedNames.add(regionName);
		}

		PrintWriter writer;
		try {
			writer = new PrintWriter(output, "ASCII");
			writer.print("Sample");
			writer.print('\t');
			writeHeader(writer, "Somatic", orderedNames);
			writer.print('\t');
			writeHeader(writer, "Tumor", orderedNames);
			writer.print('\t');
			writeHeader(writer, "Normal", orderedNames);
			writer.print('\t');
			writer.println();
	
			for(String sampleName : somaticMutRateTarget.getCountPerTypePerSample().keySet()) {				
				writer.print(sampleName);
				writeCounts(writer, sampleName, somaticMutRateTarget, orderedNames, regionSizes);
				writeCounts(writer, sampleName, tumorMutRateTarget, orderedNames, regionSizes);
				writeCounts(writer, sampleName, normalMutRateTarget, orderedNames, regionSizes);
				writer.println();
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		writer.close();
	}

	private void writeCounts(PrintWriter writer, String sampleName, MutationRateTarget mutRateTarget, List<String> orderedNames, TObjectLongHashMap<String> regionSizes) {
		TObjectLongHashMap<String> sampleMap = mutRateTarget.getCountPerTypePerSample().get(sampleName);
		writer.print('\t');
		writer.print(sampleMap.get("total"));
		writer.print('\t');
		double genomeRate = (double)sampleMap.get("total")/(double)totalSize * 1000000.0;
		writer.print(genomeRate);
		for(String name : orderedNames) {
			writer.print('\t');
			writer.print((double)sampleMap.get(name)/(double)regionSizes.get(name) * 1000000.0);
		}
		for(String name : orderedNames) {
			writer.print('\t');
			writer.print(((double)sampleMap.get(name)/(double)regionSizes.get(name) * 1000000.0)/genomeRate);
		}

		long totalCount = 0;
		for(String name : MutationRateTarget.substitutionOrder) {
			writer.print('\t');
			totalCount += mutRateTarget.getSubstitutionCountPerSample().get(sampleName).get(name);
			writer.print(mutRateTarget.getSubstitutionCountPerSample().get(sampleName).get(name));
		}
		if(totalCount != sampleMap.get("total")) {
			System.err.println("Wrong total: " + sampleName);
		}

		for(String name : MutationRateTarget.substitutionOrder) {
			writer.print('\t');
			writer.print((double)mutRateTarget.getSubstitutionCountPerSample().get(sampleName).get(name)/totalCount);
		}
		
		for(String name : MutationRateTarget.substitutionOrder) {
			for(String regionName : orderedNames) {
				writer.print('\t');
				writer.print(mutRateTarget.getSubstitutionCountPerTypePerSample().get(sampleName).get(name+'_'+regionName));
			}
		}
	}

	private void writeHeader(PrintWriter writer, String prefix, List<String> orderedNames) {
		writer.print(prefix+"Total");
		writer.print('\t');
		writer.print(prefix+"GenomeRate");
		for(String name : orderedNames) {
			writer.print("\t"+prefix);
			writer.print(name);
		}
		for(String name : orderedNames) {
			writer.print("\t"+prefix);
			writer.print(name+"_normalized");
		}
		for(String name : MutationRateTarget.substitutionOrder) {
			writer.print("\t"+prefix+'_');
			writer.print(name);
		}
		for(String name : MutationRateTarget.substitutionOrder) {
			writer.print("\t"+prefix+'_');
			writer.print(name+"_pct");
		}
		for(String name : MutationRateTarget.substitutionOrder) {
			for(String regionName : orderedNames) {
				writer.print("\t"+prefix+'_');
				writer.print(name+'_'+regionName);
			}
		}
	}
	
	public Set<Marker> getHits(Marker marker) {
		boolean hitChromo = false;
		Set<Marker> hits = new HashSet<Marker>();

		Markers intersects = predictor.intersects(marker);
		if (intersects.size() > 0) {
			for (Marker markerInt : intersects) {
				if(!(markerInt instanceof Chromosome) && !(markerInt instanceof Gene))
					hits.add(markerInt);

				if (markerInt instanceof Chromosome) {
					hitChromo = true; // OK (we have to hit a chromosome, otherwise it's an error
				} else if (markerInt instanceof Gene) {
					// Analyze Genes
					Gene gene = (Gene) markerInt;
					if(gene.numChilds() > 0) {
						// For all transcripts...
						for (Transcript tr : gene) {
	
							if (tr.intersects(marker)) { // Does it intersect this transcript?
								//hits.add(tr);
	
								for (Utr utr : tr.getUtrs())
									if (utr.intersects(marker)) hits.add(utr);
	
								for (Cds cds : tr.getCds())
									if (cds.intersects(marker)) hits.add(cds);
								//for (Exon ex : tr)
								//	if (ex.intersects(marker)) hits.add(ex);
	
								for (Intron intron : tr.introns())
									if (intron.intersects(marker)) hits.add(intron);
							}
						}
					}
					else {
						throw new RuntimeException("All genes should have transcripts");
					}
				}
			}
		}

		if (!hitChromo) throw new RuntimeException("ERROR: Out of chromosome range. " + marker);
		return hits;
	}

	public Region parseDNAse(File dnasebed) {
		Region dnase = new Region("dnase");
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(dnasebed), "ASCII"), 1*1024*1024);
			Genome genome = new Genome();
			while(true) {
				String line = reader.readLine();
				if(line == null)
					break;

				String values[] = line.split("\t");
				Marker region = new Marker(genome.getOrCreateChromosome(values[0]), Integer.parseInt(values[1]), Integer.parseInt(values[2])-1, 1, values[3]);
				dnase.add(region);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// Oh well
					e.printStackTrace();
				}
			}
		}
		
		return dnase;
	}
}
