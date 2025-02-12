# **PubMed Article Search**  

This project is a PubMed search and retrieval tool that allows users to filter articles by author and topic, and export the results to a CSV file. The tool is designed for researchers and data analysts who require quick and precise access to PubMed articles with user-defined filters.  

---

## **Project Structure**  

```plaintext
pubmed-article-search/  
├── scripts/  
│   └── get_papers_list.py  # Main script for fetching and exporting PubMed articles  
├── pyproject.toml          # Poetry configuration file for dependency management and packaging  
├── README.md               # Project documentation  
└── .gitignore              # Files and directories ignored by Git  
```

### **Explanation:**  
- **`scripts/get_papers_list.py`**: Core script that fetches PubMed articles, filters them based on user input, and exports the results to a CSV file.  
- **`pyproject.toml`**: Defines the project’s metadata, dependencies, and entry points using Poetry.  
- **`README.md`**: Provides installation instructions, usage guidelines, and an overview of the project.  
- **`.gitignore`**: Lists files and directories that should not be tracked by Git.  

---

## **Installation and Setup**  

Follow these steps to install and execute the program:  

### **1. Clone the Repository**  
   ```bash
   git clone https://github.com/Saikatsinha007/PubMed-Article-Search.git
   cd PubMed-Article-Search
   ```

### **2. Install Poetry**  
   Poetry is used for dependency management. Install it with the following commands:  
   - **Windows**:  
     ```bash
     (Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | python -
     ```
   - **macOS/Linux**:  
     ```bash
     curl -sSL https://install.python-poetry.org | python3 -
     ```

   Confirm that Poetry is installed successfully:  
   ```bash
   poetry --version
   ```

### **3. Install Project Dependencies**  
   Install all dependencies listed in `pyproject.toml` with the following command:  
   ```bash
   poetry install
   ```

### **4. Execute the Program**  
   Run the program using Poetry:  
   ```bash
   poetry run get-papers-list
   ```

---

## **Key Features**  

- **Author and Topic Filtering**: Search for PubMed articles by specific authors and topics.  
- **Date Range**: Limit the search results to articles published within a specified date range.  
- **Data Export**: Export search results as a CSV file for further analysis.  
- **User Interface**: Interactive input form built using Streamlit.  

---

## **Dependencies**  

The following tools and libraries are used in the project:  

- **[Pandas](https://pandas.pydata.org/)**: Used for data manipulation and exporting search results to CSV.  
- **[BioPython](https://biopython.org/)**: Accesses the PubMed database via NCBI’s Entrez module.  
- **[Streamlit](https://streamlit.io/)**: Provides a simple and interactive web interface for user inputs and results display.  
- **[Poetry](https://python-poetry.org/)**: Manages project dependencies and packaging, ensuring a reproducible environment.  

---

## **Contribution**  

Contributions are welcome! If you’d like to contribute to this project, please fork the repository, make your changes, and submit a pull request.  

For any issues, feel free to open a new issue in the GitHub repository.  

---

## **License**  

This project is licensed under the MIT License. See the `LICENSE` file for more details.  

---  
