package ca.mcgill.genome;

public class SequenceLengthException extends Exception {
  private static final long serialVersionUID = 7613414803307269735L;

  public SequenceLengthException() {
  }

  public SequenceLengthException(String message) {
    super(message);
  }

  public SequenceLengthException(Throwable cause) {
    super(cause);
  }

  public SequenceLengthException(String message, Throwable cause) {
    super(message, cause);
  }

  public SequenceLengthException(String message, Throwable cause, boolean enableSuppression, boolean writableStackTrace) {
    super(message, cause, enableSuppression, writableStackTrace);
  }

}
