package ca.mcgill.genome.mps.core.parser;

public class Sequence {
  private final String header;
  private final String sequence;
  private final byte[] baseQualities;

  public Sequence(String header, String sequence, byte[] baseQualities) {
    this.header = header;
    this.sequence = sequence;
    this.baseQualities = baseQualities;
  }

  public String getHeader() {
    return header;
  }

  public String getSequence() {
    return sequence;
  }

  public byte[] getBaseQualities() {
    return baseQualities;
  }

}
