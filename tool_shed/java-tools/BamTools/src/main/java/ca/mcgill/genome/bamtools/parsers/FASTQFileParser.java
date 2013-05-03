package ca.mcgill.genome.bamtools.parsers;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.nio.charset.Charset;

import org.apache.commons.io.input.CountingInputStream;
import org.obiba.bitwise.client.LargeGZIPInputStream;

public class FASTQFileParser extends SequenceFileParser {
  private static final Charset ASCII = Charset.forName("ASCII");
  private static final int BUFFER_SIZE = 128 * 1024;
  private final byte qualityOffset;

  private final CountingInputStream counter;
  private final LineNumberReader sequenceReader;

  public FASTQFileParser(File fastq, byte qualityOffset) {
    super();
    this.qualityOffset = qualityOffset;
    try {
      counter = new CountingInputStream(new BufferedInputStream(new FileInputStream(fastq), BUFFER_SIZE));
      InputStream readFrom = counter;
      if (fastq.getName().toLowerCase().endsWith(".gz")) {
        readFrom = new LargeGZIPInputStream(counter);
      }
      sequenceReader = new LineNumberReader(new InputStreamReader(readFrom, ASCII));
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * This method is not precise. It is at most a BUFFER_SIZE off.
   * @return Number of bytes read
   */
  @Override
  public long getByteCount() {
    return counter.getByteCount();
  }

  @Override
  public Sequence nextSequence() throws IOException {
    String header = sequenceReader.readLine();
    if (header == null)
      return null;
    // Remove first character
    header = header.substring(1);

    String sequence = sequenceReader.readLine();
    // Quality header
    sequenceReader.readLine();
    byte encodedQual[] = sequenceReader.readLine().getBytes(ASCII);
    byte baseQualities[] = new byte[encodedQual.length];
    int qualitySum = 0;
    for (int i = 0; i < baseQualities.length; i++) {
      baseQualities[i] = (byte) (encodedQual[i] - qualityOffset);
      qualitySum += baseQualities[i];
    }

    return new Sequence(header, sequence, baseQualities, (byte)(qualitySum/baseQualities.length));
  }

  @Override
  public void close() throws IOException {
    sequenceReader.close();
  }

}
