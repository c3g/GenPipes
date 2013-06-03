package ca.mcgill.genome.bamtools;

public abstract class DefaultTool {
	public abstract String getCmdName();
	public String getCmdUsage() {
		return "";
	}
	
	public abstract int run(String args[]);
	
	protected final void printUsageHeader() {
		System.out.println("Usage: Bamtools "+getCmdName());
	}
}
