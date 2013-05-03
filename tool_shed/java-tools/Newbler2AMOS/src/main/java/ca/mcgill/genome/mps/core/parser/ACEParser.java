/**
 * 
 */
package ca.mcgill.genome.mps.core.parser;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ca.mcgill.genome.mps.core.parser.ace.Contig;
import ca.mcgill.genome.mps.core.parser.ace.Read;
import ca.mcgill.genome.mps.core.parser.ace.Read.Direction;

/**
 *
 */
public class ACEParser {
  private static final Logger log = LoggerFactory.getLogger(ACEParser.class);
  private static Pattern DS_PATTERN = Pattern.compile("^DS +CHROMAT_FILE: +(.+) PHD_FILE: +(.+) TIME: (.+) CHEM: 454 TEMPLATE: +(.+) DIRECTION: +(.+)");

  private long nbReads;
  private long nbContigs;
  private List<Contig> contigs;

  public ACEParser(File aceFile) {
    processFileAce(aceFile);
  }

  public long getNbReads() {
    return nbReads;
  }

  public long getNbContigs() {
    return nbContigs;
  }

  public List<Contig> getContigs() {
    return contigs;
  }

  private static List<String> tokenize(String line) {
    StringTokenizer st = new StringTokenizer(line);
    List<String> tokens = new ArrayList<String>();
    while (st.hasMoreElements()) {
      String t = st.nextToken();
      tokens.add(t);
    }
    return tokens;
  }

  /**
   * Returns the position of '*' character in sequence.
   * @param sequence
   * @return
   */
  public static List<Integer> getGapsPosition(String sequence) {
    StringBuffer seqBuf = new StringBuffer(sequence);
    List<Integer> retVal = new ArrayList<Integer>();
    for (int i = 0; i < seqBuf.length(); i++) {
      char c = seqBuf.charAt(i);
      if (c == '*') {
        retVal.add(i);
      }
    }
    return retVal;
  }

  private void processFileAce(File aceFile) {
    contigs = new ArrayList<Contig>();
    try {
      BufferedReader in = new BufferedReader(new FileReader(aceFile));

      boolean done = false;

      boolean COseq = false;
      boolean CObq = false;
      boolean rdSeq = false;
      StringBuilder seqStr = new StringBuilder();
      ByteArrayOutputStream baseQualities = new ByteArrayOutputStream();
      Read rd = null;

      Contig contig = null;

      while (!done) {
        String line = in.readLine();

        if (line == null) {
          done = true;
        } else if (line.startsWith("AS")) {
          List<String> tokens = tokenize(line);
          nbContigs = Long.parseLong(tokens.get(1));
          nbReads = Long.parseLong(tokens.get(2));
        } else if (line.startsWith("CO")) {
          List<String> tokens = tokenize(line);
          contig = new Contig(tokens.get(1));
          contig.setNumBases(Integer.parseInt(tokens.get(2)));
          contig.setNumReads(Integer.parseInt(tokens.get(3)));
          contig.setNumBaseSegs(Integer.parseInt(tokens.get(4)));
          contigs.add(contig);
          COseq = true;
          CObq = false;
          rdSeq = false;
        } else if (line.startsWith("BQ")) {
          if (COseq) {
            contig.setUnTrimmedSequence(seqStr.toString());
            seqStr.setLength(0);
            COseq = false;
          }
          CObq = true;
        } else if (line.startsWith("AF")) {
          if (CObq) {
            contig.setUnTrimmedBaseQualities(baseQualities.toByteArray());
            baseQualities.reset();
            CObq = false;
          }
          List<String> tokens = tokenize(line);
          Read read = new Read(tokens.get(1));
          read.setOffset(Integer.parseInt(tokens.get(3)) - 1); // one based
          read.setComplemented(tokens.get(2).equals("C"));
          contig.addRead(read);
        } else if (line.startsWith("BS")) {
          COseq = false;
          CObq = false;
          contig.getBsLines().add(line);
        } else if (line.startsWith("RD")) {
          if (rdSeq) {
            rd.setSequence(seqStr.toString());
            seqStr.setLength(0);
            rd = null;
          }

          List<String> tokens = tokenize(line);
          rd = contig.getReads().get(tokens.get(1));
          rd.setNbPaddedBases(Integer.parseInt(tokens.get(2)));
          rd.setNbWholeReadInfoItems(Integer.parseInt(tokens.get(3)));
          rd.setNbReadTags(Integer.parseInt(tokens.get(4)));
          rdSeq = true;
        } else if (line.startsWith("QA")) {
          List<String> tokens = tokenize(line);
          rd.setQualClipStart(Integer.parseInt(tokens.get(1)));
          rd.setQualClipEnd(Integer.parseInt(tokens.get(2)));
          rd.setAlignClipStart(Integer.parseInt(tokens.get(3)));
          rd.setAlignClipEnd(Integer.parseInt(tokens.get(4)));
        } else if (line.startsWith("DS")) {
          if (rdSeq) {
            rd.setSequence(seqStr.toString());
            seqStr.setLength(0);
            rdSeq = false;
          }

          Matcher matcher = DS_PATTERN.matcher(line);
          if (matcher.matches()) {
            if (matcher.group(5).equalsIgnoreCase("fwd")) {
              // if (rd.isComplemented())
              // log.error("AF marked complemented, but direction fwd: {}", rd.getRealName());
              rd.setDirection(Direction.FWD);
            } else if (matcher.group(5).equalsIgnoreCase("rev")) {
              // if (rd.isComplemented())
              // log.error("AF marked not complemented, but direction rev: {}", rd.getRealName());
              rd.setDirection(Direction.REV);
            } else {
              log.error("Unknown direction: {}", matcher.group(5));
            }
          }
        } else if (line.length() > 0) {
          if (COseq || rdSeq) {
            seqStr.append(line);
          } else if (CObq) {
            List<String> qualities = tokenize(line);
            for (String quality : qualities) {
              baseQualities.write(Byte.parseByte(quality));
            }
          } else {
            log.error(new StringBuilder().append("extra:").append(line).toString());
          }
        }
      }
      in.close();
    } catch (Exception e) {
      log.error(e.getMessage(), e);
    }
  }
}
