# In downloads.py
from csv import DictWriter
from io import BytesIO, StringIO
from zipfile import ZipFile

from streamlit import download_button, success


class Downloads:
  """
  Handles download logic.
  """

  @staticmethod
  def csv_str(results: list[dict]) -> str:
    buffer = StringIO()
    writer = DictWriter(buffer, fieldnames=["Propriedade", "Valor"])
    writer.writeheader()
    writer.writerows(results)
    return buffer.getvalue()

  @staticmethod
  def fasta_str(title: str, sequence: str) -> str:
    return f">{title}\n{sequence}"

  @staticmethod
  def zip_download(file_name: str, csv_results: list, protein_name: str,
                   protein_sequence: str) -> None:
    # Convert data to strings
    csv_data = Downloads.csv_str(csv_results)
    fasta_data = Downloads.fasta_str(protein_name, protein_sequence)

    # Create ZIP in memory
    zip_buffer = BytesIO()
    with ZipFile(zip_buffer, "w") as zip_file:
      zip_file.writestr(f"{file_name}.csv", csv_data)
      zip_file.writestr(f"{protein_name}.fasta", fasta_data)
    zip_buffer.seek(0)

    download_button(
      label="📦Baixar Resultados (.zip)",
      data=zip_buffer,
      file_name=f"{file_name}.zip",
      mime="application/zip",
      type="primary")
