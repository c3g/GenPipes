package ca.mcgill.genome.bamtools.commands;

import gnu.trove.map.hash.TObjectIntHashMap;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.parsers.GeneMapper;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect;
import ca.mcgill.mcb.pcingola.vcf.VcfEffect.FormatVersion;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;

public class GenePerPatient extends DefaultTool {
	private int minDepth = 10;
	private GeneMapper mapper = null;
	private boolean onlyDamaging = false;

	@Override
	public String getCmdName() {
		return "genesampleimpact";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--vcf           VCFs to load");
		System.out.println("\t--mindepth      Minimum depth on which to use a Variant. (default: " + minDepth + ")");
		System.out.println("\t--mapper        Genes to keep file. Use to filter for refseq");
		System.out.println("\t--onlyDamaging  Use only damaging snps. yes,no");
	}

	@Override
	public int run(String[] args) {
		try {
			List<File> vcfs = new ArrayList<File>();
	
			for (int idx = 1; idx < args.length; idx++) {
				if (args[idx].equals("--vcf")) {
					idx++;
					vcfs.add(new File(args[idx]));
				}
				else if (args[idx].equals("--mindepth")) {
					idx++;
					minDepth = Integer.parseInt(args[idx]);
				}
				else if (args[idx].equals("--mapper")) {
					idx++;
					mapper = new GeneMapper(new File(args[idx]));
				}
				else if (args[idx].equals("--onlyDamaging")) {
					idx++;
					onlyDamaging = Boolean.parseBoolean(args[idx]);
				}
			}
	
			if (vcfs.size() == 0) {
				printUsage("vcf not set");
				return 1;
			}
			
			getImpacts(vcfs);
		} catch (Exception e) {
 			throw new RuntimeException(e);
		}

		return 0;
	}

	public Map<String, Map<String, Boolean>> testVCF(VcfFileIterator vcfParser) {
		Map<String, Map<String, Boolean>> gene2Type = new HashMap<String, Map<String, Boolean>>();

		for (VcfEntry vcfEntry : vcfParser) {
			List<VcfEffect> allEffects = vcfEntry.parseEffects(FormatVersion.FORMAT_SNPEFF_2);
			List<VcfEffect> effects = new ArrayList<VcfEffect>();

			for(VcfEffect effect : allEffects) {
				if(mapper.getTranscriptId(effect.getTranscriptId()) != null) {
					if(onlyDamaging) {
						if(vcfEntry.getInfo("dbnsfpPolyphen2_HVAR_pred").indexOf("D") != -1 || vcfEntry.getInfo("dbnsfpPolyphen2_HVAR_pred").indexOf("P") != -1 || vcfEntry.getInfoFloat("dbnsfpSIFT_score") <= 0.05f) {
							effects.add(effect);
						}
					}
					else {
						effects.add(effect);
					}
				}
			}
			
			for(VcfEffect effect : effects) {
				if(!gene2Type.containsKey(effect.getGene())) {
					gene2Type.put(effect.getGene(), new HashMap<String, Boolean>());
				}
				gene2Type.get(effect.getGene()).put(effect.getEffect().toString(), Boolean.TRUE);
			}
		}
		
		return gene2Type;
	}

	public void getImpacts(List<File> vcfs){
		Map<String, TObjectIntHashMap<String>> gene2NbTypes = new HashMap<String, TObjectIntHashMap<String>>();
		TObjectIntHashMap<String> gene2NbPatients = new TObjectIntHashMap<String>();
		Set<String> types = new HashSet<String>();
		for(File vcf : vcfs) {
			VcfFileIterator vcfParser = new VcfFileIterator(vcf.toString());
			Map<String, Map<String, Boolean>> gene2Type = testVCF(vcfParser);
			for(String gene : gene2Type.keySet()) {
				gene2NbPatients.adjustOrPutValue(gene, 1, 1);
				if(!gene2NbTypes.containsKey(gene)) {
					gene2NbTypes.put(gene, new TObjectIntHashMap<String>());
				}
				Map<String, Boolean> impacts = gene2Type.get(gene);
				for(String impact : impacts.keySet()) {
					types.add(impact);
					gene2NbTypes.get(gene).adjustOrPutValue(impact, 1, 1);
				}
			}
		}

		System.out.print("GeneName,TotalPatients");
		for(String type : types) {
			System.out.print(',');
			System.out.print(type);			
		}
		System.out.print("\n");
		
		for(String gene : gene2NbTypes.keySet()) {
			System.out.print(gene);
			System.out.print(',');
			System.out.print(gene2NbPatients.get(gene));
			for(String type : types) {
				System.out.print(',');
				System.out.print(gene2NbTypes.get(gene).get(type));			
			}
			System.out.print("\n");
		}
	}
}
