package ca.mcgill.genome.bamtools.commands;

import java.awt.Color;
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

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.labels.StandardCategoryItemLabelGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.renderer.category.BoxAndWhiskerRenderer;
import org.jfree.chart.renderer.category.StandardBarPainter;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.statistics.DefaultBoxAndWhiskerCategoryDataset;

import ca.mcgill.genome.bamtools.DefaultTool;
import ca.mcgill.genome.bamtools.mutations.CovariateAxis;

public class PlotSubstitutionRatios extends DefaultTool {
	private File ratios;
	private String outputPrefix;
	private List<String> normNames = new ArrayList<String>();
	private List<String> substitutionPctNames = new ArrayList<String>();
	private List<String> substitutionCountNames = new ArrayList<String>();
	private Map<String, List<Double>> dataNorm = new HashMap<String, List<Double>>();
	private Map<String, Map<String, Double>> dataSubstPCT = new HashMap<String, Map<String, Double>>();
	private Map<String, Map<String, Long>> dataSubstCounts = new HashMap<String, Map<String, Long>>();
	private Map<String, Map<String, Integer>> covariates = new HashMap<String, Map<String, Integer>>();
	
	@Override
	public String getCmdName() {
		return "plotsubstitutions";
	}

	@Override
	public String getCmdUsage() {
		return super.getCmdUsage();
	}

	private void printUsage(String errMsg) {
		System.out.println(errMsg);
		printUsageHeader();
		System.out.println("\t--ratios  Input ratios. Generated from the tool 'mutationrates'");
		System.out.println("\t--output  Output prefix");
		
	}

	@Override
	public int run(String[] args) {
	
		for (int idx = 1; idx < args.length; idx++) {
			if (args[idx].equals("--ratios")) {
				idx++;
				ratios = new File(args[idx]);
			}
			else if (args[idx].equals("--output")) {
				idx++;
				outputPrefix = args[idx];
			}

		}

		if(ratios == null) {
			printUsage("Missing ratios");
			return 1;
		}

		plotSubstitutionRatios();
		return 0;
	}

	public void plotSubstitutionRatios() {
		parseData();
		plotSubstitutions();
		plotMpMB();
		
	}

	private void plotMpMB() {
		DefaultBoxAndWhiskerCategoryDataset boxandwhiskercategorydataset = new DefaultBoxAndWhiskerCategoryDataset();
		for(String category : dataNorm.keySet()) {
			boxandwhiskercategorydataset.add(dataNorm.get(category), "", category);
		}

		JFreeChart jFreeChart = ChartFactory.createBoxAndWhiskerChart("Feature Mutation Rates", "Features", "Relative rate", boxandwhiskercategorydataset, true);
		jFreeChart.setBackgroundPaint(Color.WHITE);
		CategoryPlot plot = (CategoryPlot)jFreeChart.getPlot();
		plot.setBackgroundPaint(Color.WHITE);
		BoxAndWhiskerRenderer renderer =  (BoxAndWhiskerRenderer)plot.getRenderer();
		renderer.setMeanVisible(false);
//		plot.setDomainGridlinesVisible(true);
//		plot.setRangePannable(true);
//		NumberAxis numberaxis = (NumberAxis) plot.getRangeAxis();
//		numberaxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
		plot.getDomainAxis().setCategoryLabelPositions(CategoryLabelPositions.UP_45);

		try {
			ChartUtilities.saveChartAsPNG(new File(outputPrefix+"rate.png"), jFreeChart, 640, 480);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private void plotSubstitutions() {
		DefaultCategoryDataset pctCategorydataset = new DefaultCategoryDataset();
		DefaultCategoryDataset countsCategorydataset = new DefaultCategoryDataset();
		List<String> orderedSamples = new ArrayList<String>(dataSubstPCT.keySet());
		Collections.sort(orderedSamples, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				return dataSubstPCT.get(o1).get("T>A").compareTo(dataSubstPCT.get(o2).get("T>A"));
			}
		});
		for(String sample : orderedSamples) {
			for(String pctName : dataSubstPCT.get(sample).keySet()) {
				pctCategorydataset.addValue(dataSubstPCT.get(sample).get(pctName), pctName, sample);
				countsCategorydataset.addValue(dataSubstCounts.get(sample).get(pctName), pctName, sample);
			}
		}

		ChartFactory.setChartTheme(StandardChartTheme.createJFreeTheme());
		JFreeChart pctChart =  ChartFactory.createStackedBarChart("Percent of Mutation Fraction", "Samples", "Mutation Fraction", pctCategorydataset, PlotOrientation.VERTICAL, true, false, false);
		JFreeChart countsChart =  ChartFactory.createStackedBarChart("Count of Mutation Fraction", "Samples", "Mutation Fraction", countsCategorydataset, PlotOrientation.VERTICAL, true, false, false);
		pctChart.setBackgroundPaint(Color.WHITE);
		countsChart.setBackgroundPaint(Color.WHITE);

		CategoryPlot pctPlot = (CategoryPlot)pctChart.getPlot();
		CategoryPlot countsPlot = (CategoryPlot)countsChart.getPlot();
		pctPlot.setBackgroundPaint(Color.WHITE);
		countsPlot.setBackgroundPaint(Color.WHITE);

		BarRenderer renderer = (BarRenderer) pctPlot.getRenderer();
	    renderer.setShadowVisible(false);
	    renderer.setBarPainter(new StandardBarPainter());
	    renderer.setBaseItemLabelGenerator(new StandardCategoryItemLabelGenerator());
	    renderer = (BarRenderer) countsPlot.getRenderer();
	    renderer.setShadowVisible(false);
	    renderer.setBarPainter(new StandardBarPainter());
	    renderer.setBaseItemLabelGenerator(new StandardCategoryItemLabelGenerator());
	    renderer = null;
	    
	    CovariateAxis pctCovariateAxis = new CovariateAxis(covariates);
	    pctPlot.setDomainAxis(pctCovariateAxis);
	    countsPlot.setDomainAxis(pctCovariateAxis);
	    pctCovariateAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_90);
	    pctCovariateAxis.addLegends(pctChart);
	    CovariateAxis countsCovariateAxis = new CovariateAxis(covariates);
	    pctPlot.setDomainAxis(countsCovariateAxis);
	    countsPlot.setDomainAxis(countsCovariateAxis);
	    countsCovariateAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_90);
	    countsCovariateAxis.addLegends(countsChart);
		
		pctPlot.getRangeAxis().setRange(0, 1);

		try {
			ChartUtilities.saveChartAsPNG(new File(outputPrefix+".substitutionPct.png"), pctChart, 1280, 960);
			ChartUtilities.saveChartAsPNG(new File(outputPrefix+".substitutionCounts.png"), countsChart, 1280, 960);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private void parseData() {
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new InputStreamReader(new FileInputStream(ratios), "ASCII"));
			String line = reader.readLine();
			String values[] = line.split("\t");
			Map<String, Integer> normIdxs = new HashMap<String,Integer>();
			Map<String, Integer> substitutionPctIdxs = new HashMap<String,Integer>();
			Map<String, Integer> substitutionCountIdxs = new HashMap<String,Integer>();
			for(int idx=0; idx < values.length; idx++) {
				if(values[idx].endsWith("_norm") && !values[idx].startsWith("all")) {
					normNames.add(values[idx]);
					normIdxs.put(values[idx], idx);
					dataNorm.put(values[idx], new ArrayList<Double>());
				}
				else if(values[idx].endsWith("_pct") && !values[idx].startsWith("all")) {
					String nameToUse = values[idx].substring(0, values[idx].length()-"_pct".length());
					substitutionPctNames.add(nameToUse);
					substitutionPctIdxs.put(nameToUse, idx);
				}
				else if(!values[idx].endsWith("_pct") && !values[idx].startsWith("all") && values[idx].indexOf('>') != -1) {
					String nameToUse = values[idx];
					substitutionCountNames.add(nameToUse);
					substitutionCountIdxs.put(nameToUse, idx);
				}
			}

			while(true) {
				line = reader.readLine();
				if(line == null)
					break;

				if(line.startsWith("#"))
					continue;

				values = line.split("\t");
				dataSubstPCT.put(values[0], new HashMap<String, Double>());
				dataSubstCounts.put(values[0], new HashMap<String, Long>());
				covariates.put(values[0], new HashMap<String, Integer>());
//				covariates.get(values[0]).put("sex", Integer.parseInt(values[values.length-3]));
//				covariates.get(values[0]).put("age", Integer.parseInt(values[values.length-2]));
//				covariates.get(values[0]).put("grade", Integer.parseInt(values[values.length-1]));
				covariates.get(values[0]).put("sex", 0);
				covariates.get(values[0]).put("age", 0);
				covariates.get(values[0]).put("grade", 0);

				for(String normName : normNames) {
					dataNorm.get(normName).add(new Double(values[normIdxs.get(normName)]));
				}
				for(String pctName : substitutionPctNames) {
					dataSubstPCT.get(values[0]).put(pctName, new Double(values[substitutionPctIdxs.get(pctName)]));
				}
				for(String pctName : substitutionCountNames) {
					dataSubstCounts.get(values[0]).put(pctName, new Long(values[substitutionCountIdxs.get(pctName)]));
				}
			}
			reader.close();
//			reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File("D:/travail/covariate_CK.tsv")), "ASCII"));
//			while(true) {
//				line = reader.readLine();
//				if(line == null)
//					break;
//
//				values = line.split("\t");
//				covariates.put(values[0], new HashMap<String, Integer>());
//				covariates.get(values[0]).put("sex", Integer.parseInt(values[1]));
//				covariates.get(values[0]).put("age", Integer.parseInt(values[2]));
//				covariates.get(values[0]).put("grade", Integer.parseInt(values[3]));
//			}
		}catch(IOException e) {
			throw new RuntimeException(e);
		}
		finally {
			if(reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
	}
	
}
