from Query import Query

# NOT CURRENTLY IN USE

# all queries is a dict with structure:
#   {
#       databaseA:
#           {
#               query1:[QueryObj, QueryObj, QueryObj],
#               query2:[QueryObj]
#           }
#       databaseB:
#           {
#               query1:[QueryObj, QueryObj, QueryObj],
#               query2:[QueryObj]
#           }
#   }   etc.

class AllQueries:
    def __init__(self):
        all_queries = {
            "A188v1": {},
            "B73v5": {},
            "W22v2": {}
        }

    def __addHit__(self, q, genome):
        if q.query not in self.all_queries[genome]:
            self.all_queries[genome] = q.query
        self.all_queries[genome][q.query].append(q)
