---
title: "Stability analysis"
format: html
editor: visual
---

## Stability analysis

[Step1]{.underline}

Create a dummy set representing the first few rows of the original dataset.

```{r}
 
df <- data.frame( 
  cellnames = letters[1:7], 
  seed1 = c(4,0,4,0,0,4,0),
 seed2 = c(3,0,3,0,0,3,0),
 seed3 = c(14,0,14,0,0,14,0),
 seed4 = c(4,0,4,0,0,4,0),
 seed5 = c(4,5,4,2,5,4,2),
 seed6 = c(4,0,4,0,0,4,0)) 

df
```

[Step2]{.underline}

Want to change the layout of the dataframe in a wide format using pivot_long

```{r}

#cols - the columns you want to reshape to a long format (each column will change to a new entry in a row) 

#starts_with - select the columns that start with "seed" should be reshaped

#names_to - the name of the new column that will store the original columns (seed1:seed6)

#values_to - name of the new column that will store the values from the columns that are reshaped (seed1:seed6)

df_long <- df %>% 
    pivot_longer(cols = starts_with("seed"), names_to = "seed", values_to = "cluster") 

df_long
```

[Step3]{.underline}

Get all possible cell combinations by using full_join. When joining based on seed number (you are joining duplicate keys), you will get all possible combinations.

```{r}

#The full_join fuction combines two data.frames based on a common column,including both rows from both data.frames. 

#Joining df_long with df_long 

#by - a character vector specifying the columns by which the data.frames should be joined.

 comparison_df <- df_long %>% 
    full_join(df_long, by = "seed") %>% 
    filter(cellnames.x != cellnames.y) 
  

```

[Step4]{.underline}

Removing the same observations (eg cell a vs b == cell b vs cell a). To do this you need to sort the dataframe so that every observation will be in the same order. Then see which rows are duplicated and then remove these rows.

```{r}
 
#apply() function: is used to apply a specified function to either rows or columns of an array, matrix, or data frame.
 
#In this case, apply(comparison_df, 1, sort) applies the sort() function to each row of the data frame df.
 
#Sorting Rows:By specifying MARGIN = 1, we indicate that we want to apply the function to rows (i.e., across each row).

#The sort() function sorts the elements within each row in ascending order.

 #Result:The resulting object sorted will be a matrix where each row contains the sorted values from the corresponding row in the original data frame.
 
#When you apply t() to an object (such as a matrix or a data frame), it flips the rows and columns.
 
#Essentially, it swaps the rows with the columns, creating a new object with the dimensions reversed.
 
 #now that each row is arranged in the same order you can see which rows are duplicates (eg  cell ab vs cell ba)
 
 sorted <- apply(comparison_df,1,sort)
 transposed <- t(sorted) 
 
 
comparison_df <- comparison_df[!duplicated(transposed), ]
  

```

[Step5]{.underline}

If cell.x and cell.y are in the same cluster for the same seed, a 1 will be allocated in the 'same_cluster' column, else a 0

```{r}

  comparison_df <- comparison_df %>% 
    mutate(same_cluster = ifelse(cluster.x == cluster.y, 1,0)) 
```

[Step6]{.underline}

Selecting the cellnames, seed_value and same_cluster and changing the layout(changing seeds from observations to variables).This is so that the same_cluster values for each cell pair in a specific seed can be added across the rows.

```{r}

comparison_matrix <- comparison_df %>% 
    select(cellnames.x,cellnames.y, seed, same_cluster) %>% 
    pivot_wider(names_from = "seed", values_from = "same_cluster") 

```

[Step7]{.underline}

Summing the amount of times that each cell pair cluster together for each seed number.

```{r}

trial <- comparison_matrix %>% 
  mutate(total = rowSums(comparison_matrix[,3:8])) %>%
  select(cellnames.x,cellnames.y,total) 

#calculating the percentage for cell pair clustering 

 trial <- trial %>% 
   mutate(percentage = (total*100)/6)


```
