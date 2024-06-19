import Layout from '../components/Layout';

export default function Example() {
  // Render example data page with information and links
  return (
    <Layout>
      <div>
        <h1 className="h3 mb-3 text-gray-800">Example Data Page</h1>
        <h2 className="h4 mb-3 text-gray-800">Bacteroides thetaiotaomicron - plasmid</h2>
        <p>
          Bacteroides thetaiotaomicron is a Gram-negative bacterium that is known for its ability to degrade complex polysaccharides in the human gut.
          The following small dataset includes one B. thetaiotaomicron strain, focusing only on the plasmid of the orgnaims.
          This is especially usefull for testing the interface.
          The data was computed based on the study by <a href="https://www.nature.com/articles/s41467-020-17348-5">Ryan et al. (2020) <i class="fas fa-external-link-alt"></i></a>.
          <br />The dataset is available for download here: <a href="/example/bacteroidesplasmid/dataset.zip"><i className="fas fa-download"></i></a>
          <br />The corresponding report from TSS-Captur is available here: <a href="/example/bacteroidesplasmid/report/Interface/overview.html">[Click here]</a>
        </p>
        <h2 className="h4 mb-3 text-gray-800">Streptomyces coelicolor</h2>

        <p>
          Streptomyces coelicolor is a soil-dwelling bacterium that is known for its complex life cycle and its ability to produce a wide range of secondary metabolites. The following dataset includes one S. coelicolor strain <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=100226">A3(2) <i class="fas fa-external-link-alt"></i></a>.
          The data was originally from the study by  <a href="https://www.nature.com/articles/ncomms11605">Jeong et al. (2016) <i class="fas fa-external-link-alt"></i></a>.
          <br />The dataset is available for download here: <a href="/example/scoelicolor/dataset.zip"><i className="fas fa-download"></i></a>
          <br />The corresponding report from TSS-Captur is available here: <a href="/example/scoelicolor/report/Interface/overview.html">[Click here]</a>
        </p>
        <h2 className="h4 mb-3 text-gray-800">Campylobacter jejuni</h2>

        <p>
          Campylobacter jejuni, a Gram-negative bacterium, is known to be one of the leading causes of bacterial gastroenteritis in humans. Its adaptability to a wide range of hosts and conditions, as well as its comparatively small genome of about 1.6 Mb, make it a good choice for the study of prokaryotic transcriptional mechanisms.
          The following dataset includes four C. jejuni strains: <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=192222">NCTC11168 <i class="fas fa-external-link-alt"></i></a>, <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=354242">81-176 <i class="fas fa-external-link-alt"></i></a>, <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=407148">81116 <i class="fas fa-external-link-alt"></i></a> and <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=195099">RM1221 <i class="fas fa-external-link-alt"></i></a>. These strains were subjected to a high resolution transcriptome analysis by <a href="https://doi.org/10.1371/journal.pgen.1003495">Dugar et al. (2013) <i class="fas fa-external-link-alt"></i></a>.
          <br />The dataset is available for download here: <a href="/example/cjejuni/dataset.zip"><i className="fas fa-download"></i></a>
          <br />The corresponding report from TSS-Captur is available here: <a href="/example/cjejuni/report/Interface/overview.html">[Click here]</a>
        </p>
      </div>
    </Layout>
  );
}