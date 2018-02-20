package fr.ens.biologie.genomique.eoulsan.modules.preprocessing;

import fr.ens.biologie.genomique.eoulsan.EoulsanException;
import fr.ens.biologie.genomique.eoulsan.EoulsanRuntimeDebug;
import fr.ens.biologie.genomique.eoulsan.bio.GFFEntry;
import fr.ens.biologie.genomique.eoulsan.bio.io.GFFReader;

import org.junit.*;
import org.junit.rules.TemporaryFolder;
import org.junit.rules.ExpectedException;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class TestParser {

    private final File annotations =
        new File("src/test/files/testAnnotation.gff");
    private final File genesRef = new File("src/test/files/genesTable.tsv");
    private final File transRef = new File("src/test/files/transTable.tsv");
    private final File exonsRef = new File("src/test/files/exonsTable.tsv");
    private final File extransRef = new File("src/test/files/extransTable.tsv");

    @Before public void setUp() throws Exception {
        EoulsanRuntimeDebug.initDebugEoulsanRuntime();
    }

    @Rule public TemporaryFolder folder = new TemporaryFolder();

    @Rule public ExpectedException thrown = ExpectedException.none();

    @Test public void testGetMetadataGene() {
        try {

            GFFReader annotationReader = new GFFReader(annotations);
            SCFeatureMetadata gene = new SCFeatureMetadata();

            for (GFFEntry anno : annotationReader) {
                String type = anno.getType();

                if (type.equals("gene") && gene.getId().equals("")) {
                    gene = FeaturesMetadataExtractorModule
                        .getMetadata(anno, false, "ID", "MT", "ERCC");
                }

                if ( ! gene.getId().equals("")) {
                    break;
                }

            }

            Assert.assertEquals("6", gene.getChromosome());
            Assert.assertEquals("ENSMUSG00000000058", gene.getId());
            Assert.assertEquals(17281185, gene.getStart());
            Assert.assertEquals(17289115, gene.getEnd());

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    @Test public void testGetMetadataTranscript() {
        try {

            GFFReader annotationReader = new GFFReader(annotations);
            SCFeatureMetadata transcript = new SCFeatureMetadata();


            for (GFFEntry anno : annotationReader) {
                String type = anno.getType();


                if (type.equals("transcript") && transcript.getId().equals("")) {
                    transcript = FeaturesMetadataExtractorModule
                        .getMetadata(anno, false, "ID", "MT", "ERCC");
                }

                if ( ! transcript.getId().equals("")) {
                    break;
                }

            }


            Assert.assertEquals("6", transcript.getChromosome());
            Assert.assertEquals("ENSMUST00000000058", transcript.getId());
            Assert.assertEquals(17281185, transcript.getStart());
            Assert.assertEquals(17289115, transcript.getEnd());

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    @Test public void testGetMetadataExon() {
        try {

            GFFReader annotationReader = new GFFReader(annotations);
            SCFeatureMetadata exon = new SCFeatureMetadata();

            for (GFFEntry anno : annotationReader) {
                String type = anno.getType();


                if (type.equals("exon") && exon.getId().equals("")) {
                    exon = FeaturesMetadataExtractorModule
                        .getMetadata(anno, false, "Name", "MT", "ERCC");
                }

                if(! exon.getId().equals("")) {
                    break;
                }

            }

            Assert.assertEquals("6", exon.getChromosome());
            Assert.assertEquals("ENSMUSE00000713825", exon.getId());
            Assert.assertEquals(17281185, exon.getStart());
            Assert.assertEquals(17281510, exon.getEnd());
            Assert.assertEquals("ENSMUST00000000058", exon.getTranscript());


        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    @Test public void testSaveMetadata() {

        try {

            GFFReader annotationReader = new GFFReader(annotations);

            Map<String, SCFeatureMetadata> features = new HashMap();

            for (GFFEntry anno : annotationReader) {

                SCFeatureMetadata feature = FeaturesMetadataExtractorModule
                    .getMetadata(anno, false,
                        anno.getType().equals("exon") ? "Name" : "ID", "MT",
                        "ERCC");

                FeaturesMetadataExtractorModule.saveFeature(features, feature);
            }

            Assert.assertEquals(21, features.size());

            // Reinitialize
            features = new HashMap();
            annotationReader = new GFFReader(annotations);

            for (GFFEntry anno : annotationReader) {

                SCFeatureMetadata feature = FeaturesMetadataExtractorModule
                    .getMetadata(anno, false,
                        anno.getType().equals("exon") ? "Parent" : "ID", "MT",
                        "ERCC");

                FeaturesMetadataExtractorModule.saveFeature(features, feature);
            }

            Assert.assertEquals(10, features.size());
            Assert.assertEquals(17281185, features.get("ENSMUST00000000058").getStart());
            Assert.assertEquals(17289115, features.get("ENSMUST00000000058").getEnd());
            Assert.assertEquals(10664, features.get("ENSMUST00000000058").getLength());

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Test public void testExtractMetadata() {

        try {

            File outputGenes = folder.newFile();
            File outputTrans = folder.newFile();
            File outputExons = folder.newFile();
            File outputExtrans = folder.newFile();

            FeaturesMetadataExtractorModule
                .extractMetadata(annotations, "MT", "ERCC", "gene", "ID", false,
                    outputGenes);
            FeaturesMetadataExtractorModule
                .extractMetadata(annotations, "MT", "ERCC", "transcript", "ID",
                    false, outputTrans);
            FeaturesMetadataExtractorModule
                .extractMetadata(annotations, "MT", "ERCC", "exon", "Name",
                    false, outputExons);
            FeaturesMetadataExtractorModule
                .extractMetadata(annotations, "MT", "ERCC", "exon", "Parent",
                    false, outputExtrans);

            // Genes
            BufferedReader refReader =
                new BufferedReader(new FileReader(genesRef));
            BufferedReader resReader =
                new BufferedReader(new FileReader(outputGenes));
            Assert.assertEquals(resReader.readLine(), refReader.readLine());

            // Transcripts
            refReader = new BufferedReader(new FileReader(transRef));
            resReader = new BufferedReader(new FileReader(outputTrans));
            Assert.assertEquals(resReader.readLine(), refReader.readLine());

            // Exons/Transcripts
            refReader = new BufferedReader(new FileReader(extransRef));
            resReader = new BufferedReader(new FileReader(outputExtrans));
            Assert.assertEquals(resReader.readLine(), refReader.readLine());

            // Exons
            refReader = new BufferedReader(new FileReader(exonsRef));
            resReader = new BufferedReader(new FileReader(outputExons));
            Assert.assertEquals(resReader.readLine(), refReader.readLine());

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    @Test public void testExtractionError() throws Exception {

        //test type
        thrown.expect(EoulsanException.class);

        //test message
        thrown.expectMessage("Feature exon does not contain a ID attribute");

        File output = folder.newFile();
        FeaturesMetadataExtractorModule
            .extractMetadata(annotations, "MT", "ERCC", "exon", "ID", false,
                output);

    }
}

