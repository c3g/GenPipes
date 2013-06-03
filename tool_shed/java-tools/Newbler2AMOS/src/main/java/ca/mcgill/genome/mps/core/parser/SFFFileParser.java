package ca.mcgill.genome.mps.core.parser;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.charset.Charset;

import org.apache.commons.io.input.CountingInputStream;

/**
 * Taken from the roche manual File Formats GS-FLX-System-Software-Manual-SW-Manual-Overview-File-Formats
 */
public class SFFFileParser {
  private final Charset ASCII = Charset.forName("ASCII");
  private String version;
  private int nbOfReads;
  private int readIdx;
  // keep in case we find a use for the indexes
  @SuppressWarnings("unused")
  private long indexOffset;
  // keep in case we find a use for the indexes
  @SuppressWarnings("unused")
  private int indexLength;
  private short numberOfFlowsPerRead;
  private int flowgramBytesPerFlow;
  private CountingInputStream counter;
  private DataInputStream sffInput;

  public SFFFileParser(File sffFile) {
    super();
    try {
      counter = new CountingInputStream(new BufferedInputStream(new FileInputStream(sffFile)));
      sffInput = new DataInputStream(counter);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    setup();
  }

  private void setup() {
    try {
      readIdx = 0;

      int magicNumber = sffInput.readInt();
      if (magicNumber != 0x2E736666) {
        throw new RuntimeException("Magic Number failed: " + magicNumber);
      }

      byte rawVersion[] = new byte[4];
      int read = sffInput.read(rawVersion);
      if (read != rawVersion.length) {
        throw new RuntimeException("Didn't read the whole version");
      }

      version = new String(rawVersion);
      indexOffset = sffInput.readLong();
      indexLength = sffInput.readInt();
      nbOfReads = sffInput.readInt();
      short headerLength = sffInput.readShort();
      short keyLength = sffInput.readShort();
      numberOfFlowsPerRead = sffInput.readShort();
      byte flowgramFormatCode = sffInput.readByte();
      if (flowgramFormatCode != 1) {
        throw new RuntimeException("Unknown flowgram code: " + flowgramFormatCode);
      }
      flowgramBytesPerFlow = 2; // Because flowgramFormatCode == 1
      // uint16_t, where the floating point flowgram value is encoded as "(int) round(value * 100.0)", and decoded as
      // "(storedvalue * 1.0 / 100.0)".

      char flowChars[] = new char[numberOfFlowsPerRead];// [numberOfFlowsPerRead]
      for (int i = 0; i < numberOfFlowsPerRead; i++) {
        flowChars[i] = (char) sffInput.readByte();
      }

      char keySequence[] = new char[keyLength]; // [key_length]
      for (int i = 0; i < keyLength; i++) {
        keySequence[i] = (char) sffInput.readByte();
      }

      int withoutPadding = 31 + numberOfFlowsPerRead + keyLength;
      for (int i = 0; i < (headerLength - withoutPadding); i++) {
        if (sffInput.readByte() != 0) {
          throw new RuntimeException("Padding is wrong. wo: " + withoutPadding + " with: " + headerLength);
        }
      }
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  public long getByteCount() {
    return counter.getByteCount();
  }

  public Sequence nextSequence() throws IOException {
    if (readIdx == nbOfReads)
      return null;

    short readHeaderLength = sffInput.readShort();
    short nameLength = sffInput.readShort();
    int numberOfBases = sffInput.readInt();
    short clipQualLeft = sffInput.readShort();
    short clipQualRight = sffInput.readShort();
    short clipAdapterLeft = sffInput.readShort();
    short clipAdapterRight = sffInput.readShort();

    byte rawName[] = new byte[nameLength];
    int read = sffInput.read(rawName);
    if (read != rawName.length) {
      throw new RuntimeException("Didn't read the read name");
    }
    String name = new String(rawName, ASCII);

    int withoutPadding = 16 + nameLength;
    for (int i = 0; i < (readHeaderLength - withoutPadding); i++) {
      if (sffInput.readByte() != 0) {
        throw new RuntimeException("Read Header padding is wrong. wo: " + withoutPadding + " with: " + readHeaderLength);
      }
    }

    short flowgramValues[] = new short[numberOfFlowsPerRead];
    for (int i = 0; i < numberOfFlowsPerRead; i++) {
      flowgramValues[i] = sffInput.readShort();
    }

    byte flowIndexPerBase[] = new byte[numberOfBases];
    read = sffInput.read(flowIndexPerBase);
    if (read != flowIndexPerBase.length) {
      throw new RuntimeException("Didn't read the whole flowIndexPerBase");
    }

    byte rawBases[] = new byte[numberOfBases];
    read = sffInput.read(rawBases);
    if (read != rawBases.length) {
      throw new RuntimeException("Didn't read the whole bases");
    }
    String bases = new String(rawBases, ASCII);

    byte qualityScores[] = new byte[numberOfBases];
    read = sffInput.read(qualityScores);
    if (read != qualityScores.length) {
      throw new RuntimeException("Didn't read the whole qualityScores");
    }

    int length = numberOfFlowsPerRead * flowgramBytesPerFlow + 3 * numberOfBases;
    int totalLength = length;
    if (length % 8 != 0) {
      totalLength = length - (length % 8) + 8;
    }

    for (int i = 0; i < (totalLength - length); i++) {
      if (sffInput.readByte() != 0) {
        throw new RuntimeException("Read padding is wrong. wo: " + length + " with: " + totalLength);
      }
    }

    readIdx++;
    return new SoftClippedSequence(name, bases, qualityScores, clipQualLeft, clipQualRight, clipAdapterLeft, clipAdapterRight);
  }

  public void close() throws IOException {
    try {
      sffInput.close();
    } finally {
      sffInput = null;
      counter = null;
    }
  }

  public String getVersion() {
    return version;
  }

  public int getNbOfReads() {
    return nbOfReads;
  }

  public int getFlowgramBytesPerFlow() {
    return flowgramBytesPerFlow;
  }

}
