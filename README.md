# 🌟 PubMed Article Search

🔍 **PubMed Article Search** is a tool that allows users to search and retrieve PubMed articles, filter by author and topic, and export the results to a CSV file. It is ideal for researchers and data analysts who need precise and filtered access to PubMed data.

---

## 📁 Project Structure

```plaintext
pubmedsearch/
├── pubmed_module.py       # Core module with functions to search PubMed
├── streamlit_pubmed.py    # Streamlit web app for interactive article search
├── cli_pubmed.py          # Command-line interface for fetching and exporting articles
├── pyproject.toml         # Poetry configuration file for dependency management and packaging
├── README.md              # Project documentation
├── LICENSE                # Project license (MIT)
└── .env                   # Environment variables (e.g., ENTREZ_API_KEY)
```

### ✨ Explanation
- **`pubmed_module.py`**: Contains the core logic to search PubMed articles using the NCBI Entrez API and returns results as a pandas DataFrame.
- **`streamlit_pubmed.py`**: Provides an interactive web interface using Streamlit for searching PubMed articles and downloading results as CSV.
- **`cli_pubmed.py`**: Implements a command-line interface using `argparse` to search PubMed articles and export the results to a CSV file. This script is also registered as an executable command.
- **`pyproject.toml`**: Defines the project’s metadata, dependencies, and entry points using Poetry. It registers the CLI command as `get-papers-list`.
- **`README.md`**: Contains installation instructions, usage guidelines, and an overview of the project.
- **`LICENSE`**: The project is licensed under the MIT License.
- **`.env`**: Used to store sensitive environment variables such as the NCBI Entrez API key.

---

## 🛠️ Installation and Setup

Follow these steps to install and execute the program:

### 1. Clone the Repository
```bash
git clone https://github.com/Saikatsinha007/PubMed-Article-Search.git
cd PubMed-Article-Search/pubmedsearch
```

### 2. Install Poetry

Poetry is used for dependency management. Install it using one of the following commands:

- **Windows**:
  ```bash
  (Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | python -
  ```
- **macOS/Linux**:
  ```bash
  curl -sSL https://install.python-poetry.org | python3 -
  ```

✅ Verify Poetry is installed:
```bash
poetry --version
```

### 3. Install Project Dependencies

Install all dependencies defined in `pyproject.toml` by running:
```bash
poetry install
```

### 4. Configure Environment Variables

Create a `.env` file in the project root (if not already present) and add your NCBI Entrez API key:
```ini
ENTREZ_API_KEY=your_entrez_api_key_here
```

---

## ▶️ Execution

### Command-Line Interface (CLI)

After installing dependencies, you can run the CLI tool using Poetry:
```bash
poetry run get-papers-list --start_date 2012-03-01 --end_date 2022-12-31 --authors "Smith, Doe" --topics "Cancer, Genomics"
```
This command will fetch the PubMed articles based on your input filters and export the results to a CSV file.

### Streamlit Web App

To run the interactive web interface, execute:
```bash
poetry run streamlit run streamlit_pubmed.py
```
This will open the Streamlit app in your browser, where you can perform searches and download the CSV results directly.

---

## 🌟 Key Features

- **Author and Topic Filtering**: Search for PubMed articles by specific authors and topics.
- **Date Range**: Limit search results to articles published within a specified date range.
- **Data Export**: Export search results as a CSV file for further analysis.
- **Multiple Interfaces**: Choose between the command-line tool for quick searches or the Streamlit web app for an interactive experience.

---

## 📚 Dependencies

This project uses the following tools and libraries:

- **[Pandas](https://pandas.pydata.org/)**: For data manipulation and exporting search results to CSV.
- **[BioPython](https://biopython.org/)**: Provides access to PubMed via NCBI’s Entrez module.
- **[Streamlit](https://streamlit.io/)**: Offers a simple, interactive web interface for input and result display.
- **[Poetry](https://python-poetry.org/)**: Manages project dependencies and packaging.

---

## 🤝 Contribution

Contributions are welcome! 🎉 To contribute:

1. Fork the repository.
2. Create a new branch with your changes.
3. Submit a pull request with a detailed explanation of your changes.

For any issues, please open a new issue in the GitHub repository.

---

## 📝 License

This project is licensed under the MIT License. See the `LICENSE` file for details.
```

