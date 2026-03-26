# VRP Network Solver - Case Study

Este repositório contém uma solução personalizada para o **Problema de Roteirização de Veículos (VRP)** aplicada a uma rede viária real. O projeto utiliza uma biblioteca própria desenvolvida em Python, que integra o poder de otimização do **Google OR-Tools** com a manipulação de dados espaciais.

O exemplo contido no notebook utiliza a rede de transporte de **Lisboa** para demonstrar a criação de rotas otimizadas entre diferentes pontos (nós) da cidade.

## 🚀 Funcionalidades

* **Modelagem de Rede:** Conversão de arquivos de nós (`nodes`) e conexões (`links`) em um grafo direcionado usando `NetworkX`.
* **Cálculo de Matriz de Distância:** Implementação da fórmula de Haversine para precisão em coordenadas geográficas.
* **Otimização com OR-Tools:** Resolução de problemas de logística buscando a minimização de distância e tempo.
* **Integração Geoespacial:** Exportação dos resultados em `GeoDataFrame` (GeoPandas) para fácil visualização e análise de mapas.
* **Visualização:** Suporte para plotagem das rotas geradas sobre a malha urbana.

## 📂 Estrutura do Projeto

* `vrpnetwork.py`: A biblioteca principal. Contém a classe `NetworkVRP` que gerencia a rede, calcula as matrizes e invoca o solver do OR-Tools.
* `VRP_LIS.ipynb`: Notebook demonstrativo que carrega os dados de Lisboa, configura os parâmetros de roteirização e exibe os resultados.
* `nodes_lis2.csv` / `links_lis2.csv`: Base de dados da rede viária de Lisboa (necessários para rodar o exemplo).
* `roteirizacao.jpg`: Exemplo visual do output das rotas otimizadas.

## 🛠️ Pré-requisitos

Para rodar este projeto, você precisará das seguintes bibliotecas:

```bash
pip install ortools networkx pandas geopandas matplotlib shapely
```

## 📖 Como Usar

1.  **Instancie a rede:**
    ```python
    from vrpnetwork import NetworkVRP
    network = NetworkVRP()
    network.load_data("nodes_lis2.csv", "links_lis2.csv")
    ```

2.  **Defina os pontos de demanda:**
    Selecione os IDs dos nós que representam o depósito e os clientes na rede.

3.  **Resolva o problema:**
    O método interno utiliza o OR-Tools para encontrar a sequência lógica que minimiza o custo total da operação.

4.  **Visualize os resultados:**
    O notebook demonstra como converter a solução em geometrias `LineString` e `Point` para visualização em mapas.

## 📊 Exemplo de Resultado

Abaixo, uma representação das rotas calculadas sobre a malha de Lisboa:

![Roteirização em Lisboa](roteirizacao.jpeg)

## 👤 Autor

**Fábio Emanuel de Souza Morais**
*Mestre em Engenharia de Transportes (USP) | Engenheiro Aeronáutico (ITA)*

