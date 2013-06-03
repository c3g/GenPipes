package ca.mcgill.genome.parser;

public class Sequence {
  private final String header;
  private final String sequence;
  private final byte[] baseQualities;
  private final byte avgQuality;

  public Sequence(String header, String sequence, byte[] baseQualities, byte avgQuality) {
    this.header = header;
    this.sequence = sequence;
    this.baseQualities = baseQualities;
    this.avgQuality = avgQuality;
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

  public byte getAvgQuality() {
    return avgQuality;
  }

}
