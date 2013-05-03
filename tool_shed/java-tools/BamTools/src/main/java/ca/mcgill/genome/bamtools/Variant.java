package ca.mcgill.genome.bamtools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ca.mcgill.mcb.pcingola.util.Gpr;

public class Variant {
	private final String chromosome;
	private final int position;
	private final String gene;
	private String refAllele;
	private final List<String> altAlleles;
	private final List<SampleGenotype> sampleGenotypes;
	private final Map<String, SampleGenotype> sampleToGenotype;
	
	public Variant(String chromosome, int position, String refAllele, String altAlleles) {
		this(chromosome, position, refAllele, altAlleles, "");
	}

	public Variant(String chromosome, int position, String refAllele, String altAlleles, String gene) {
		this.chromosome = chromosome;
		this.position = position;
		this.gene = gene;
		this.refAllele = refAllele;
		String alleles[] = altAlleles.split(",");
		this.altAlleles = Arrays.asList(alleles);
		this.sampleGenotypes = new ArrayList<SampleGenotype>();
		this.sampleToGenotype = new HashMap<String, SampleGenotype>(); 
	}

	
	public static Variant read(String line) {
		List<String> values = new ArrayList<String>(Arrays.asList(Gpr.split(line, '\t')));
		Variant variant = new Variant(values.remove(0), Gpr.parseIntSafe(values.remove(0)), values.remove(0), values.remove(0), values.remove(0));
		while(values.size() > 3) {
			variant.addSampleGenotype(values.remove(0), values.remove(0), values.remove(0), Double.parseDouble(values.remove(0)), Integer.parseInt(values.remove(0)));
		}
		
		return variant;
	}

	public String write() {
		StringBuilder sb = new StringBuilder();
		sb.append(chromosome).append('\t');
		sb.append(position).append('\t');
		sb.append(gene).append('\t');
		sb.append(refAllele).append('\t');
		for(int idx=0; idx < altAlleles.size();idx++) {
			if(idx > 0)
				sb.append(',');
			sb.append(altAlleles.get(idx));
		}
		sb.append('\t');
		
		for(SampleGenotype sg : sampleGenotypes) {
			sg.write(sb);
		}
		
		return sb.toString();
	}

	public void simplifyIndels() {
		String newRef = refAllele;
		String newAlleles[] = new String[altAlleles.size()];
		for(int idx=0; idx < altAlleles.size(); idx++) {
			newAlleles[idx] = altAlleles.get(idx);
		}

		while(newRef.length() > 1) {
			boolean isOk = true;
			for(String newAllele: newAlleles) {
				if(newAllele.length() == 1) {
					isOk = false;
					break;
				}
				if(!isOk)
					break;
			}
			
			char lastNuc = refAllele.charAt(refAllele.length()-1);
			for(String newAllele: newAlleles) {
				if(newAllele.charAt(newAllele.length()-1) != lastNuc) {
					isOk = false;
					break;
				}
			}
			
			if(!isOk) {
				break;
			}

			newRef = newRef.substring(0, newRef.length()-1);
			for(int idx=0; idx < altAlleles.size(); idx++) {
				newAlleles[idx] = newAlleles[idx].substring(0, newAlleles[idx].length()-1);
			}
		}

		refAllele = newRef;
		for(int idx=0; idx < newAlleles.length; idx++) {
                        altAlleles.set(idx, newAlleles[idx]);
                }
	}

	public boolean isIndel() {
		if(refAllele.length() > 1) {
			return true;
		}
		
		for(String altAllele : altAlleles) {
			if(altAllele.length() > 1) {
				return true;
			}
		}
		return false;
	}

	public String getChromosome() {
		return chromosome;
	}

	public int getPosition() {
		return position;
	}

	public String getGene() {
		return gene;
	}

	public String getRefAllele() {
		return refAllele;
	}

	public List<String> getAltAlleles() {
		return altAlleles;
	}
	
	public Map<String, SampleGenotype> getSampleToGenotype() {
		return sampleToGenotype;
	}

	public List<SampleGenotype> getSampleGenotypes() {
		return sampleGenotypes;
	}

	public void addSampleGenotype(String sampleName, String alleleA, String alleleB, double qual, int clr) {
		SampleGenotype sg = new SampleGenotype(sampleName, alleleA, alleleB, qual, clr);
		sampleGenotypes.add(sg);
		if(sampleToGenotype.put(sampleName, sg) != null) {
			throw new RuntimeException("Sample added twice: "+ sampleName);
		}
	}

	public static class SampleGenotype {
		private final String sampleName;
		private final String alleleA;
		private final String alleleB;
		private final double qual;
		private final int clr;

		public SampleGenotype(String sampleName, String alleleA, String alleleB) {
			this(sampleName, alleleA, alleleB, -1, -1);
		}

		public SampleGenotype(String sampleName, String alleleA, String alleleB, double qual, int clr) {
			super();
			this.sampleName = sampleName;
			this.alleleA = alleleA;
			this.alleleB = alleleB;
			this.qual = qual;
			this.clr = clr;
		}

		private void write(StringBuilder sb) {
			sb.append(sampleName).append('\t');
			sb.append(alleleA).append('\t');
			sb.append(alleleB).append('\t');
			sb.append(qual).append('\t');
			sb.append(clr).append('\t');
		}

		public boolean isVariant(String refAllele) {
			return !refAllele.equals(alleleA) || !isHom();
		}

		public boolean isHom() {
			return alleleA.equals(alleleB);
		}

		public String getSampleName() {
			return sampleName;
		}

		public String getAlleleA() {
			return alleleA;
		}

		public String getAlleleB() {
			return alleleB;
		}

		public double getQual() {
			return qual;
		}

		public int getClr() {
			return clr;
		}
		
	}
}
