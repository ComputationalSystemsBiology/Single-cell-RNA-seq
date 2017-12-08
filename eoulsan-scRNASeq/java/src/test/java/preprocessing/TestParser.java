package preprocessing;


import fr.ens.biologie.genomique.eoulsan.EoulsanException;
import fr.ens.biologie.genomique.eoulsan.Globals;
import fr.ens.biologie.genomique.eoulsan.annotations.LocalOnly;
import fr.ens.biologie.genomique.eoulsan.bio.GFFEntry;
import fr.ens.biologie.genomique.eoulsan.bio.io.GFFReader;
import fr.ens.biologie.genomique.eoulsan.bio.io.GTFReader;
import fr.ens.biologie.genomique.eoulsan.core.*;
import fr.ens.biologie.genomique.eoulsan.data.DataFile;
import org.junit.Test;


import java.io.*;
import java.util.HashMap;
import java.util.Map;


public class TestParser {

    /**
     * Extract and Write features metadata from annotation file
     *
     * @param annotations the annotation file
     * @param mtTag       tag for mitochondrial features
     * @param spikeTag    tag for spike in features
     * @param genomicType genomic type of the features to consider
     * @param attributeId Id to design the feature
     * @param gtfFormat   boolean indicating if the annotation file is gtf or not
     * @param outFile     output file to write data
     * @throws IOException      if encouters problem with input or output file
     * @throws EoulsanException if features as no Id
     */
    private void extractMetadata(File annotations, final String mtTag, final String spikeTag,
                                 final String genomicType, final String attributeId, final boolean gtfFormat,
                                 File outFile) throws IOException, EoulsanException {

        try (final GFFReader annotationReader =
                     gtfFormat ? new GTFReader(new FileInputStream(annotations)) : new GFFReader(new FileInputStream(annotations));
             BufferedWriter out = new BufferedWriter(new FileWriter(outFile))) {

            // Initialise gathering variables
            Map<String, SCFeatureMetadata> features = new HashMap<>();
            Map<String, String> parents = new HashMap<>();

            //Read annotation file
            for (final GFFEntry anno : annotationReader) {

                // Treat non genic case, saving gene ID
                if (!genomicType.equals("gene")) {
                    parents = saveGeneId(anno, parents, gtfFormat);
                }

                // Extract information
                if (genomicType.equals(anno.getType())) {

                    String type;

                    // get Id
                    final String featureId = anno.getAttributeValue(attributeId);

                    if (featureId == null) {
                        throw new EoulsanException("Feature " + genomicType + " does not contain a "
                                + attributeId + " attribute");
                    }

                    // get length
                    final int start = anno.getStart();
                    final int end = anno.getEnd();
                    final int length = end - start;

                    // get feature type

                    if (anno.getSeqId().equals(mtTag)) {
                        type = "mitochondrial";
                    } else if (anno.getSeqId().equals(spikeTag)) {
                        type = "spike";
                    } else {
                        type = "other";
                    }

                    // Test if object is referenced (sparse features) else create it
                    if (!features.containsKey(featureId)) {
                        SCFeatureMetadata feature = new SCFeatureMetadata(genomicType);
                        feature.setLength(length);
                        feature.setType(type);
                        feature.setId(featureId);

                        features.put(featureId, feature);
                    } else {
                        SCFeatureMetadata feature = features.get(featureId);
                        final int oldLength = feature.getLength();
                        feature.setLength(oldLength + length);
                        features.put(featureId, feature);
                    }
                }
            }
            for (String id : features.keySet()) {
                SCFeatureMetadata feature = features.get(id);
                if (!feature.isGene()) {
                    feature.setGene(parents.get(id));
                }
            }

            // Write headers in output file
            out.write( "Type" + "\t" + "Length");
            if (!genomicType.equals("gene")) {
                out.write("\t" + "GeneID");
            }
            // Write metadata in output file
            out.newLine();
            for (String id : features.keySet()) {
                features.get(id).printMetadata(out);
                out.newLine();
            }
        }
    }

    /**
     * Save genes id in a already existing Map
     *
     * @param anno      a GFFentry object
     * @param map       Map to save data
     * @param gtfFormat boolean indicating whether annotations are from a gtf file or not
     * @return Map completed with key : feature id, value : corresponding gene id
     */
    private Map<String, String> saveGeneId(GFFEntry anno, Map<String, String> map, boolean gtfFormat) {
        if (gtfFormat && anno.getType().matches("(.)*transcript(.)*")) {
            map.put(anno.getAttributeValue("transcript_id"), anno.getAttributeValue("gene_id"));
            return map;
        }
        if (anno.getType().matches("(.)*transcript(.)*")) {
            map.put(anno.getAttributeValue("ID"), anno.getAttributeValue("Parent"));
            return map;
        }
        return map;
    }



    @Test
    public void testParser() {
        File annotations = new File("../../files/testAnnotation.gff");

        File b = new File("../../files/genesTable.tsv");


        try {
            extractMetadata(annotations, "MT", "ERCC", "exon", "Parent", false, b);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }
}

