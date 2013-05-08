package ca.mcgill.genome.bamtools;

import java.util.LinkedHashMap;
import java.util.Map;

import ca.mcgill.genome.bamtools.commands.AlleleSpecificMethyl;
import ca.mcgill.genome.bamtools.commands.BAMHomology;
import ca.mcgill.genome.bamtools.commands.BaseFrequency;
import ca.mcgill.genome.bamtools.commands.BinnedReadCounter;
import ca.mcgill.genome.bamtools.commands.CompareFrequency;
import ca.mcgill.genome.bamtools.commands.CompareSamples;
import ca.mcgill.genome.bamtools.commands.DepthOfCoverage;
import ca.mcgill.genome.bamtools.commands.FilterDuplicates;
import ca.mcgill.genome.bamtools.commands.GenderEstimator;
import ca.mcgill.genome.bamtools.commands.GenePerPatient;
import ca.mcgill.genome.bamtools.commands.MotifCISBPRNA;
import ca.mcgill.genome.bamtools.commands.MutationRates;
import ca.mcgill.genome.bamtools.commands.PlotSubstitutionRatios;
import ca.mcgill.genome.bamtools.commands.ReplaceReadGroups;
import ca.mcgill.genome.bamtools.commands.VCF2MAF;
import ca.mcgill.genome.bamtools.commands.VariantAlleleFrequency;
import ca.mcgill.genome.bamtools.commands.VariantValidator;

public class Bamtools {
	private static Map<String, DefaultTool> tools = new LinkedHashMap<String, DefaultTool>();
	
	static {
		DefaultTool tool;
		tool = new AlleleSpecificMethyl();
		tools.put(tool.getCmdName(), tool);
		tool = new BAMHomology();
		tools.put(tool.getCmdName(), tool);
		tool = new BaseFrequency();
		tools.put(tool.getCmdName(), tool);
		tool = new CompareFrequency();
		tools.put(tool.getCmdName(), tool);
		tool = new CompareSamples();
		tools.put(tool.getCmdName(), tool);
		tool = new DepthOfCoverage();
		tools.put(tool.getCmdName(), tool);
		tool = new FilterDuplicates();
		tools.put(tool.getCmdName(), tool);
		tool = new GenderEstimator();
		tools.put(tool.getCmdName(), tool);
		tool = new GenePerPatient();
		tools.put(tool.getCmdName(), tool);
		tool = new MotifCISBPRNA();
		tools.put(tool.getCmdName(), tool);
		tool = new MutationRates();
		tools.put(tool.getCmdName(), tool);
		tool = new PlotSubstitutionRatios();
		tools.put(tool.getCmdName(), tool);
		tool = new ReplaceReadGroups();
		tools.put(tool.getCmdName(), tool);
		tool = new BinnedReadCounter();
		tools.put(tool.getCmdName(), tool);
		tool = new VariantAlleleFrequency();
		tools.put(tool.getCmdName(), tool);
		tool = new VariantValidator();
		tools.put(tool.getCmdName(), tool);
		tool = new VCF2MAF();
		tools.put(tool.getCmdName(), tool);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 1) {
			printUsageAndExit();
		}
		
		DefaultTool tool = tools.get(args[0]);
		if(tool == null) {
			printUsageAndExit();
		}

		System.exit(tool.run(args));
	}

	private static void printUsageAndExit() {
		System.out.println("Usage: Bamtools <cmd>");
		for(String toolName : tools.keySet()) {
			System.out.print('\t');
			System.out.println(toolName);
		}
		System.exit(1);
	}
}
