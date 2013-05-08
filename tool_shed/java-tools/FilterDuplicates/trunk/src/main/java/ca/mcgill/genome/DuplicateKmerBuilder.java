package ca.mcgill.genome;

import gnu.trove.TObjectIntHashMap;

import org.obiba.bitwise.BitVector;

import ca.mcgill.genome.parser.Sequence;

public class DuplicateKmerBuilder {
  private static final int ESTIMATED_MAX_ITEMS = 200000000;
  private static final int DUPLICATE_OFFSET = 1024;
  private transient TObjectIntHashMap<BitVector> nbKnownSequencesFound;
  private int kmerSize = 20;
  private int offset = 15;
  private long totalNbReads = 0;
  private boolean pairedSequence = false;

  public void processSequence(Sequence sequence) {
    if (sequence.getSequence().length() < kmerSize + offset)
      return;
    if(pairedSequence && sequence.getSequence().length() < (sequence.getSequence().length()/2)+offset+(kmerSize/2))
      return;

    if (nbKnownSequencesFound == null)
      nbKnownSequencesFound = new TObjectIntHashMap<BitVector>((int) (ESTIMATED_MAX_ITEMS / 0.75 + 1.0), 0.75f);

    totalNbReads++;
    BitVector twoBitSequence;
    try {
      twoBitSequence = computeKmer(sequence.getSequence());
    } catch (SequenceLengthException e) {
      throw new RuntimeException("Kmer too long, or offset to big: " + sequence.getHeader());
    }
    // Keep this for debugging purposes
    // if(pair1.toString().equals(pair2.toString())) {
    //   System.err.println("They are equal: " + sequence.getHeader());
    // }
    if (!nbKnownSequencesFound.containsKey(twoBitSequence)) {
      nbKnownSequencesFound.put(twoBitSequence, sequence.getAvgQuality());
    } else {
      int value = nbKnownSequencesFound.get(twoBitSequence);
      if (value >= DUPLICATE_OFFSET) {
        value = value - DUPLICATE_OFFSET;
      }

      if(value < sequence.getAvgQuality()) {
        nbKnownSequencesFound.put(twoBitSequence, DUPLICATE_OFFSET + sequence.getAvgQuality());
      }
    }
  }

  public BitVector computeKmer(String sequence) throws SequenceLengthException {
    int bitIdx = 0;
    BitVector twoBitSequence = new BitVector(kmerSize * 3);

    int nucIdx = offset;
    StringBuilder pair1 = new StringBuilder();
    StringBuilder pair2 = new StringBuilder();
    StringBuilder sb = pair1;
    for (int kmerIdx = 0; kmerIdx < kmerSize; kmerIdx++) {
      if(sequence.length() <= nucIdx) {
        throw new SequenceLengthException("Offset: " + nucIdx);
      }
      char nucleotide = sequence.charAt(nucIdx);
      sb.append(nucleotide);
      switch (nucleotide) {
      case 'A':
        // 000
        bitIdx += 3;
        break;
      case 'C':
        // 001
        bitIdx += 2;
        twoBitSequence.set(bitIdx);
        bitIdx++;
        break;
      case 'G':
        // 010
        bitIdx++;
        twoBitSequence.set(bitIdx);
        bitIdx += 2;
        break;
      case 'T':
        // 011
        bitIdx++;
        twoBitSequence.set(bitIdx);
        bitIdx++;
        twoBitSequence.set(bitIdx);
        bitIdx++;
        break;
      case 'N':
      case '-':
      case '*':
        // 100
        twoBitSequence.set(bitIdx);
        bitIdx += 3;
        break;
      default:
        throw new RuntimeException("Unknown base: " + nucleotide);
      }

      nucIdx++; 
      if(pairedSequence && nucIdx == offset+(kmerSize/2)) {
        nucIdx = (sequence.length()/2)+offset;
        sb = pair2;
      }
    }
    
    return twoBitSequence;
  }

  public boolean getKmerCountAndMark(String sequence, int quality) throws SequenceLengthException {
    BitVector twoBitSequence = computeKmer(sequence);
    int value = nbKnownSequencesFound.get(twoBitSequence);
    if(value - DUPLICATE_OFFSET <= quality) {
      nbKnownSequencesFound.put(twoBitSequence, Integer.MAX_VALUE);
      return true;
    }
    
    return false;
  }

  public TObjectIntHashMap<BitVector> getNbKnownSequencesFound() {
    return nbKnownSequencesFound;
  }

  public long getTotalNbReads() {
    return totalNbReads;
  }

  public void setTotalNbReads(long totalNbReads) {
    this.totalNbReads = totalNbReads;
  }

  public int getKmerSize() {
    return kmerSize;
  }

  public void setKmerSize(int kmerSize) {
    this.kmerSize = kmerSize;
  }

  public int getOffset() {
    return offset;
  }

  public void setOffset(int offset) {
    this.offset = offset;
  }

  public boolean isPairedSequence() {
    return pairedSequence;
  }

  public void setPairedSequence(boolean pairedSequence) {
    this.pairedSequence = pairedSequence;
  }
  
  public static class KmerCount {
    private int count;
    private byte maxQuality;
    public int getCount() {
      return count;
    }
    public void setCount(int count) {
      this.count = count;
    }
    public byte getMaxQuality() {
      return maxQuality;
    }
    public void setMaxQuality(byte maxQuality) {
      this.maxQuality = maxQuality;
    }
  }
}
