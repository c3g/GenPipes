package ca.mcgill.genome.mps.core.parser;

import java.util.Arrays;

public class SoftClippedSequence extends Sequence {
  /**
   * At what index the good sequence starts This is not a length, it's an index.
   */
  private final short clipQualLeft;
  /**
   * At what index the good sequence stops. This is not a length, it's an index. If my sequence is 500bp long and I what
   * to trim the last 490, clipQualRight will be == to 10, not 490.
   */
  private final short clipQualRight;
  private final short clipAdapterLeft;
  private final short clipAdapterRight;
  private final String unTrimmedSequence;
  private final byte[] unTrimmedBaseQualities;

  public SoftClippedSequence(String header, String unTrimmedSequence, byte[] unTrimmedBaseQualities, short clipQualLeft, short clipQualRight,
      short clipAdapterLeft, short clipAdapterRight) {
    super(header, unTrimmedSequence.substring(clipQualLeft + clipAdapterLeft - 1, clipQualRight + clipAdapterRight), Arrays.copyOfRange(unTrimmedBaseQualities,
        clipQualLeft + clipAdapterLeft - 1, clipQualRight + clipAdapterRight));

    this.unTrimmedSequence = unTrimmedSequence;
    this.unTrimmedBaseQualities = unTrimmedBaseQualities;
    this.clipQualLeft = clipQualLeft;
    this.clipQualRight = clipQualRight;
    this.clipAdapterLeft = clipAdapterLeft;
    this.clipAdapterRight = clipAdapterRight;
  }

  public short getClipQualLeft() {
    return clipQualLeft;
  }

  public short getClipQualRight() {
    return clipQualRight;
  }

  public short getClipAdapterLeft() {
    return clipAdapterLeft;
  }

  public short getClipAdapterRight() {
    return clipAdapterRight;
  }

  public String getUnTrimmedSequence() {
    return unTrimmedSequence;
  }

  public byte[] getUnTrimmedBaseQualities() {
    return unTrimmedBaseQualities;
  }

}
