//
//  filter.hpp
//  ReadGnomad
//
//  Created by Gang Peng on 4/26/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#ifndef filter_hpp
#define filter_hpp

#include <stdio.h>

#include <map>
#include <string>
#include <vector>
/*
 * -filterLarger AF 0.2 ... // filter out all variants with AF > 0.2
 * -filterSmaller AF 0.01 ...   // filter out all variants with AF < 0.01
 * -filterNotEqual variant_type snv ...  // include snv only
 * the first two for numeric variables and the last for strings
 */
int checkFilter(std::map<std::string, std::vector<std::string> > cmLine,
                 std::map<std::string, double> & filterLarger,
                 std::map<std::string, double> & filterSmaller,
                 std::map<std::string, std::vector<std::string> > & filterNotEqual);


#endif /* filter_hpp */
