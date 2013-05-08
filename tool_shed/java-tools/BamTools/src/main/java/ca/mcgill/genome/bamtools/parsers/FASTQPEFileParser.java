package ca.mcgill.genome.bamtools.parsers;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 */
public class FASTQPEFileParser extends SequenceFileParser {
  private FASTQFileParser parser1;
  private FASTQFileParser parser2;

  public FASTQPEFileParser(File fastq1, File fastq2, byte qualityOffset) {
    super();
    parser1 = new FASTQFileParser(fastq1, qualityOffset);
    parser2 = new FASTQFileParser(fastq2, qualityOffset);
  }

  @Override
  public Sequence nextSequence() throws IOException {
    Sequence sequence1 = parser1.nextSequence();
    Sequence sequence2 = parser2.nextSequence();

    if (sequence1 == null)
      return null;

    // CASAVA 1.6:
    // @HWUSI-EAS100R:6:73:941:1973#0/1
    // @HWUSI-EAS100R:6:73:941:1973#0/2

    // CASAVA 1.8:
    // @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    // @EAS139:136:FC706VJ:2:2104:15343:197393 2:Y:18:ATCACG

    String header1 = sequence1.getHeader();
    String header2 = sequence2.getHeader();

    if (header1.trim().contains(" ") && header2.trim().contains(" ")) {
      // CASAVA 1.8
      header1 = header1.substring(0, header1.indexOf(' '));
      header2 = header2.substring(0, header2.indexOf(' '));
    } else {
      // CASAVA 1.6
      header1 = header1.substring(0, sequence1.getHeader().length() - 2);
      header2 = header2.substring(0, sequence2.getHeader().length() - 2);
    }

    if (!header1.equals(header2)) {
      throw new RuntimeException("The 2 fastq files don't have the same ordering of reads.");
    }

    byte qualities[] = Arrays.copyOf(sequence1.getBaseQualities(), sequence1.getBaseQualities().length + sequence2.getBaseQualities().length);
    System.arraycopy(sequence2.getBaseQualities(), 0, qualities, sequence1.getBaseQualities().length, sequence2.getBaseQualities().length);
    
    int qualitySum = 0;
    for (int i = 0; i < qualities.length; i++) {
      qualitySum += qualities[i];
    }
    
    Sequence retVal = new Sequence(header1, sequence1.getSequence() + sequence2.getSequence(), qualities, (byte)(qualitySum/qualities.length));
    return retVal;
  }

  @Override
  public long getByteCount() {
    return parser1.getByteCount() + parser2.getByteCount();
  }

  @Override
  public void close() throws IOException {
    IOException ioException = null;
    try {
      parser1.close();
    } catch (IOException e) {
      ioException = e;
    }
    try {
      parser2.close();
    } catch (IOException e) {
      ioException = e;
    }

    if (ioException != null) {
      throw new IOException("One or both of the parsers had an exception.", ioException);
    }
  }

}
