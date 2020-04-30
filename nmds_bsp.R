library(vegan)
library(tidyverse)

#Dummy-Daten erstellen. Die Struktur sollte einer klassischen Community-Matrix entsprechen:
#          Art1  Art2  Art3  Art4...
# Probe1    2     7     7     12
# Probe2    5     12    2     2
# Probe3    8     13    2     3
# .
# .
# .
#Die kuenstliche Beispielcommunity besteht aus 5 Arten die an 60 Orten gezaehlt wurden. Der Eintrag in jeder Art steht
#fuer die Anzahl der gezaehlten Individuen, welche jeweils aus einem Zufallsgenerator gezogen wird.
community <- tibble(OrtsID = 1:60,
                    temperatur = rep(c('kalt','mittel','warm'), each = 20),
                    Art1 = sample(1:200,60),
                    Art2 = sample(12:200,60),
                    Art3 = sample(1:100,60),
                    Art4 = sample(1:300,60),
                    Art5 = sample(1:200,60))

#Fuehre nmds aus: metaMDS wird eine reine Community-Matrix uebergeben. Enthaelt der Datensatz z.B. Zusatzinfos wie ueber den Ort,
#muessen diese entfernt werden (in diesem Beispiel gibt es eine Spalte "OrtsID")
nmdsCommunity <- community %>% 
  select(-OrtsID, -temperatur) %>% 
  metaMDS(., k=2) #k = 2 steht fuer die Anzahl der Dimensionen auf die reduziert werden soll. 2 Dimensionen bedeutet, dass du
#anschliessend einen normalen Plot machen kannst, bei 3 Dimensionen muesstest du einen 3D-Plot erstellen, der oft schwer zu interpretieren ist

#Extrahiere die Scores von jedem Ort (der Score gibt die x-und y-Koordinaten fuer jeden Ort in der runterskalierten Ebene an)
#Anschliessend werden die Zusatzinformationen ueber die Orte wieder angehaengt, um Gruppen/treatments voneinander unterscheiden 
#zu koennen im Plots. Die Reihenfolge der Scores entspricht der Reihenfolge der treatments in "community"
communityScores <- as_tibble(scores(nmdsCommunity)) %>% 
  mutate(OrtsID = community$OrtsID,
         temperatur = community$temperatur)
ggplot(communityScores, aes(x = NMDS1, y = NMDS2, colour = temperatur)) +
  geom_point() +
  theme(panel.background = element_rect(fill = NA, colour = "#666666")) +
  annotate('text', x = 0.2, y = -0.24, label = paste('Stress:', round(nmdsCommunity$stress,3)))
#Der Plot zeigt logischerweise keine Cluster, da die Daten zufallsgeneriert sind. 
#Der Stress-Wert sollte immer mit dargestellt oder in der Bildunterschrift auftauchen, da er 
#ein Indikator dafuer ist wie gut die "wahren" Abstaende der Proben in der auf 2D-runterskalierten
#Darstellung erhalten sind. Der Stress von 0.22 in diesem Plot ist ziemlich schlecht

