import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import LineString,Point

import math
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

def invertecoords(lista):
    lista2 = []
    for coord in lista:
        temp = (coord[1],coord[0])
        lista2.append(temp)
    return lista2


def haversine(coord1, coord2):
    """Calcula a distância entre duas coordenadas na esfera."""
    R = 6371  # Raio da Terra em km
    lat1, lon1 = coord1
    lat2, lon2 = coord2

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    
    a = math.sin(dlat / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    return R * c  # Retorna a distância em km

class Node:
    """Classe para representar um nó com ID, latitude e longitude."""
    def __init__(self, node_id, latitude, longitude):
        self.id = node_id
        self.latitude = latitude
        self.longitude = longitude

class Link:
    """Classe para representar um link entre dois nós com atributos adicionais."""
    def __init__(self, link_id, node_a, node_b, speed, length, time, direction):
        self.id = link_id
        self.node_a = node_a
        self.node_b = node_b
        self.speed = speed
        self.length = length
        self.time = time
        self.direction = direction
        self.par = (node_a,node_b)

class Network:
    """Classe para criar e gerenciar a rede usando NetworkX."""
    def __init__(self):
        self.graph = nx.DiGraph()
        self.nodes = {}
        self.links = []

    def add_node(self, node):
        """Adiciona um nó ao grafo e ao dicionário de nós."""
        self.nodes[node.id] = node
        self.graph.add_node(node.id, pos=(node.latitude, node.longitude))

    def add_link(self, link):
        """Adiciona um link ao grafo, respeitando a direção."""
        if link.direction == 1:
            self.graph.add_edge(link.node_a, link.node_b, id=link.id, 
                                speed=link.speed, length=link.length, time=link.time)
        elif link.direction == -1:
            self.graph.add_edge(link.node_b, link.node_a, id=link.id, 
                                speed=link.speed, length=link.length, time=link.time)
        elif link.direction == 0:
            self.graph.add_edge(link.node_a, link.node_b, id=link.id, 
                                speed=link.speed, length=link.length, time=link.time)
            self.graph.add_edge(link.node_b, link.node_a, id=link.id, 
                                speed=link.speed, length=link.length, time=link.time)
        self.links.append(link)

    def load_data(self, nodes_file, links_file):
        """Carrega nós e links a partir de arquivos CSV."""
        nodes_df = (pd.read_csv(nodes_file)
                      .assign(Longitude = lambda x:x.Longitude/1000000,
                               Latitude = lambda x:x.Latitude/1000000))
        links_df = pd.read_csv(links_file)

        # Adicionar nós
        for _, row in nodes_df.iterrows():
            node = Node(row['ID'], row['Latitude'], row['Longitude']    )
            self.add_node(node)

        # Adicionar links
        for _, row in links_df.iterrows():
            link = Link(row['ID'], row['a_node'], row['b_node'], 
                        row['Vel'], row['Lenght'], row['time_ab'], row['Dir'])
            self.add_link(link)

    def shortest_path(self, source, target, weight='length'):
        """Calcula o caminho mínimo entre dois nós."""
        try:
            path = nx.shortest_path(self.graph, source=source, target=target, weight=weight)
            print(f"Caminho mínimo ({weight}): {path}")

            # Calcular a soma dos pesos (length)
            total_length = sum(
                self.graph[u][v]['length'] for u, v in zip(path[:-1], path[1:])
            )
            total_time = sum(
                self.graph[u][v]['time'] for u, v in zip(path[:-1], path[1:])
            )
            print(f"Soma total do comprimento (length) do caminho: {total_length}")

            links = []
            for i in range(len(path) - 1):
                node_a = path[i]
                node_b = path[i + 1]
                link_data = self.graph.get_edge_data(node_a, node_b)
                link_id = link_data['id']
                links.append(link_id)
                

            return path, total_length, total_time, links



         
        except nx.NetworkXNoPath:
            print(f"Não há caminho disponível de {source} para {target}.")
            return None

    def show_links_in_path(self, path):
        """Exibe os links selecionados no caminho mínimo."""
        if not path:
            print("Nenhum caminho encontrado.")
            return

        print(f"Caminho encontrado: {path}")
        for i in range(len(path) - 1):
            node_a = path[i]
            node_b = path[i + 1]
            link_data = self.graph.get_edge_data(node_a, node_b)
            link_id = link_data['id']
            print(f"Link ID={link_id} de {node_a} para {node_b}")

    def visualize_network(self):
        """Visualiza o grafo e o caminho encontrado usando Matplotlib."""
        pos = nx.get_node_attributes(self.graph, 'pos')
        plt.figure(figsize=(10, 7))
        nx.draw(self.graph, pos, with_labels=True, node_color='skyblue', node_size=700, edge_color='gray')
        plt.title("Visualização da Rede")
        plt.show()


    def visualize_shortest_path(self, path):
        """Visualiza apenas o caminho mínimo encontrado."""
        if not path:
            print("Nenhum caminho encontrado para visualizar.")
            return

        pos = nx.get_node_attributes(self.graph, 'pos')
        plt.figure(figsize=(10, 7))

        # Desenhar o grafo original em cinza
        nx.draw(self.graph, pos, with_labels=True, node_color='lightgray', node_size=500, edge_color='lightgray')

        # Desenhar o caminho mínimo em destaque
        edges_in_path = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
        nx.draw_networkx_edges(self.graph, pos, edgelist=edges_in_path, edge_color='red', width=2)
        nx.draw_networkx_nodes(self.graph, pos, nodelist=path, node_color='skyblue', node_size=700)

        plt.title("Caminho Mínimo na Rede")
        plt.show()

   


    def export_shortest_path_to_shapefile(self, path, output_file):
        """Exporta o caminho mínimo para um shapefile com uma única feature LineString."""
        if not path:
            print("Nenhum caminho encontrado para exportar.")
            return

        # Construir a geometria LineString com as coordenadas do caminho
        #for node in path:
        #    print(node)
   
        coordinates = [self.graph.nodes[node]['pos'] for node in path]
        #print(coordinates)

        invertida = invertecoords(coordinates)

        line = LineString(invertida)
        #print(line)

        # Calcular os atributos: Length, Time, node_a, node_b
        total_length = sum(
            self.graph[path[i]][path[i + 1]]['length'] for i in range(len(path) - 1)
        )
        total_time = sum(
            self.graph[path[i]][path[i + 1]]['time'] for i in range(len(path) - 1)
        )
        node_a = path[0]  # Primeiro nó
        node_b = path[-1]  # Último nó

        # Criar um GeoDataFrame com uma única feature
        gdf = gpd.GeoDataFrame(
            [{
                'node_a': node_a,
                'node_b': node_b,
                'Length': total_length,
                'Time': total_time,
                'geometry': line
            }],
            crs="EPSG:4326"
        )

        # Exportar para o arquivo shapefile
        gdf.to_file(output_file, driver='ESRI Shapefile')

        print(f"Caminho mínimo exportado para {output_file} com projeção EPSG:4326")
        


class NetworkVRP(Network):
    def create_distance_matrix(self, depot, clients):
        """Cria a matriz de distâncias entre o depósito e os clientes."""
        nodes = [depot] + clients  # Primeiro nó é o depósito
        num_nodes = len(nodes)

        # Criar a matriz de distâncias usando Haversine
        distance_matrix = [[0] * num_nodes for _ in range(num_nodes)]
        for i in range(num_nodes):
            for j in range(num_nodes):
                if i != j:
                    coord1 = self.graph.nodes[nodes[i]]['pos']
                    coord2 = self.graph.nodes[nodes[j]]['pos']
                    distance_matrix[i][j] = haversine(coord1, coord2)
        return distance_matrix

    def solve_vrp(self, depot, clients):
        """Resolve o VRP simples usando OR-Tools."""
        # Criar a matriz de distâncias
        distance_matrix = self.create_distance_matrix(depot, clients)

        # Número de nós (depósito + clientes)
        num_nodes = len(distance_matrix)

        # Criação do roteador
        manager = pywrapcp.RoutingIndexManager(num_nodes, 1, 0)  # 1 veículo, começando do nó 0 (depósito)
        routing = pywrapcp.RoutingModel(manager)

        # Função de custo (distância entre os nós)
        def distance_callback(from_index, to_index):
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return int(distance_matrix[from_node][to_node] * 1000)  # Multiplicar por 1000 para evitar arredondamento

        transit_callback_index = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        # Configurar a busca (heurísticas)
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC

        # Resolver o problema
        solution = routing.SolveWithParameters(search_parameters)

        if solution:
            return self.get_solution(manager, routing, solution)
        else:
            print("Não foi encontrada nenhuma solução.")
            return None

    def get_solution(self, manager, routing, solution):
        """Extrai a rota e a distância total da solução."""
        route = []
        index = routing.Start(0)  # Começa pelo depósito
        route_distance = 0

        while not routing.IsEnd(index):
            node = manager.IndexToNode(index)
            route.append(node)
            previous_index = index
            index = solution.Value(routing.NextVar(index))
            route_distance += routing.GetArcCostForVehicle(previous_index, index, 0)

        route.append(manager.IndexToNode(index))  # Voltar ao depósito
        print(f"Rota: {route}")
        print(f"Distância total: {route_distance / 1000} km")
        return route, route_distance / 1000  # Converter para km


    def trata_sol_VRP(self,depot,clients,route):
        path = []
        for u in route:
            if u == 0:
                path.append(depot)
            else:
                path.append(clients[u-1])
        path_final = []
        dist_final = 0
        distancias = []
        customers = []
        times = []
        links = []
        for n in range(len(path)-1):
            path2, distancia, tempo, link = self.shortest_path(path[n], path[n+1], weight='length')
            distancias.append(distancia)
            times.append(tempo)
            links.append(link)
            customers.append((path[n], path[n+1]))
            if n == 0:
                path_final = path_final + path2[:]
            else:
                path_final = path_final + path2[1:]
            dist_final = dist_final + distancia
        
        return path_final, distancias, customers, times, links, dist_final

    def caminhos(self,path):
        path_final = []
        dist_final = 0
        distancias = []
        customers = []
        times = []
        links = []
        for n in range(len(path)-1):
            path2, distancia, tempo, link = self.shortest_path(path[n], path[n+1], weight='length')
            distancias.append(distancia)
            times.append(tempo)
            links.append(link)
            customers.append((path[n], path[n+1]))
            if n == 0:
                path_final = path_final + path2[:]
            else:
                path_final = path_final + path2[1:]
            dist_final = dist_final + distancia
        
        return path_final, distancias, customers, times, links, dist_final


    def export_shortest_customers_to_shapefile(self,customers,distancias,times,links,arquivo):
        geometries = []
        dados = []
        dist_Acum = 0.0
        tempo_Acum = 0.0
        id_inicial = customers[0][0]
        dados.append({ 'sequencia': 1,
                    'id_customers': id_inicial,
                    'distancia': 0.0,
                    'tempo': 0.0,
                    'dist_Acum': dist_Acum,
                    'tempo_Acum': tempo_Acum,
                    'links' : []
                        })
        # Coordenada do ponto (longitude, latitude)
        coord = self.graph.nodes[id_inicial]['pos']
        geometries.append(Point(coord[1],coord[0]))

        for i, point in enumerate(customers):
            id = customers[i][1]
            dist_Acum = distancias[i] + dist_Acum
            tempo_Acum = times[i] + tempo_Acum
            dados.append( { 'sequencia': i + 2,
                            'id_customers': id,
                            'distancia': distancias[i],
                            'tempo': times[i],
                            'dist_Acum': dist_Acum,
                            'tempo_Acum': tempo_Acum,
                            'links': links[i]
                          })
            coord = self.graph.nodes[id]['pos']
            geometries.append(Point(coord[1],coord[0]))

         # Criar um GeoDataFrame com as geometrias e atributos
        gdf = gpd.GeoDataFrame(dados, geometry=geometries, crs="EPSG:4326")

        # Exportar o GeoDataFrame para um shapefile
        gdf.to_file(arquivo, driver='ESRI Shapefile')    