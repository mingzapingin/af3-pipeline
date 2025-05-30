<!-- PROJECT LOGO -->
<br />
<div align="center">

<h3 align="center">Alphafold3 Pipeline</h3>

  <p align="center">
    A reproducible pipeline for converting AlphaFold3 output CIFs to PDB, scoring with Prodigy and pDockQ.
    <br />
    <a href="https://github.com/mingzapingin/af3-pipeline"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/mingzapingin/af3-pipeline">View Demo</a>
    &middot;
    <a href="https://github.com/mingzapingin/af3-pipeline/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    &middot;
    <a href="https://github.com/mingzapingin/af3-pipeline/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
A reproducible pipeline for converting AlphaFold3 output CIFs to PDB, scoring with Prodigy and pDockQ.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple example steps.

### Prerequisites

You need Python 3.11+ and pip.

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/mingzapingin/af3-pipeline.git
   ```
2. Install dependencies
   ```sh
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

  '''sh
    python src/process_result.py \
      --input_folder path/to/cif_files \
      --output_folder path/to/output_folder \
      --temperature 37
  '''
    
  The script automatically uses:

  '''sh
  vendor/prodigy_prot/src/prodigy_prot/Modified_predict_IC.py
  vendor/pdockq.py
  '''

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Roadmap

- [ ] Process AF3's confidence values
    - [ ] ptm
    - [ ] iptm

See the [open issues](https://github.com/mingzapingin/af3-pipeline/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the project_license. See `LICENSE` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Kimhun Tuntikawinwong - ktuntika@bu.edu

Project Link: [https://github.com/mingzapingin/af3-pipeline](https://github.com/mingzapingin/af3-pipeline)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments
* [AlphaFold3](https://github.com/deepmind/alphafold)
* [PRODIGY](https://github.com/haddocking/prodigy)
* [pDockQ](https://github.com/fteufel/alphafold-peptide-receptors)
* [README_Template](https://github.com/othneildrew/Best-README-Template)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
