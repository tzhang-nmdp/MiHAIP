package analysis;
/**
 * Identifies which individual a variant or transcript came from.
 * Reconstructed from usage (original source not found in any tree or jar):
 * MergeTool uses Patient.Donor / Patient.Reciepor, and SNPdata/Transcript call getID().
 */
public enum Patient {
    Donor("D"),
    Reciepor("R");
    private final String id;
    Patient(String id){
        this.id = id;
    }
    public String getID(){
        return id;
    }
}
