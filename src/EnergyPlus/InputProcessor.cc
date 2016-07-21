// C++ Headers
#include <algorithm>
#include <istream>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <fstream>
#include <iostream>
#include <InputProcessor.hh>
#include <nlohmann/json.hpp>

// ObjexxFCL Headers
#include <ObjexxFCL/Backspace.hh>
#include <ObjexxFCL/Array.functions.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/gio.hh>
#include <ObjexxFCL/stream.functions.hh>
#include <ObjexxFCL/string.functions.hh>

// EnergyPlus Headers
#include <CommandLineInterface.hh>
#include <InputProcessor.hh>
#include <DataIPShortCuts.hh>
#include <DataOutputs.hh>
#include <DataPrecisionGlobals.hh>
#include <DataSizing.hh>
#include <DataStringGlobals.hh>
#include <DataSystemVariables.hh>
#include <DisplayRoutines.hh>
#include <SortAndStringUtilities.hh>

using json = nlohmann::json;

json EnergyPlus::InputProcessor::jdf = json();
json EnergyPlus::InputProcessor::schema = json();
IdfParser EnergyPlus::InputProcessor::idf_parser = IdfParser();
State EnergyPlus::InputProcessor::state = State();
std::ostream * EnergyPlus::InputProcessor::echo_stream = nullptr;

json IdfParser::decode( std::string const & idf, json const & schema ) {
	bool success = true;
	return decode( idf, schema, success );
}

json IdfParser::decode( std::string const & idf, json const & schema, bool & success ) {
	json root;
	success = true;
	if ( idf.empty() ) return root;

	size_t index = 0;
	root = parse_idf( idf, index, success, schema );
	return root;
}

std::string IdfParser::encode( json const & root, json const & schema ) {
	std::string encoded;

	for ( auto obj = root.begin(); obj != root.end(); ++obj ) {
		const auto & legacy_idd = schema[ "properties" ][ obj.key() ][ "legacy_idd" ][ "fields" ];
		for ( auto obj_in = obj.value().begin(); obj_in != obj.value().end(); ++obj_in ) {
			encoded += obj.key();
			int skipped_fields = 0;
			for ( int i = 0; i < legacy_idd.size(); i++ ) {
				std::string entry = legacy_idd[ i ];
				if ( obj_in.value().find( entry ) == obj_in.value().end() ) {
					if ( entry == "name" ) encoded += ",\n  " + obj_in.key();
					else skipped_fields++;
					continue;
				}
				for ( int j = 0; j < skipped_fields; j++ ) encoded += ",\n  ";
				skipped_fields = 0;
				encoded += ",\n  ";
				auto const & val = obj_in.value()[ entry ];
				if ( val.is_string() ) encoded += val.get < std::string >();
				else encoded += std::to_string( val.get < double >() );
			}

			if ( obj_in.value().find( "extensions" ) == obj_in.value().end() ) {
				encoded += ";\n\n";
				continue;
			}

			auto & extensions = obj_in.value()[ "extensions" ];
			for ( int extension_i = 0; extension_i < extensions.size(); extension_i++ ) {
				auto & cur_extension_obj = extensions[ extension_i ];
				auto & extensible = schema[ "properties" ][ obj.key() ][ "legacy_idd" ][ "extensibles" ];
				for ( int i = 0; i < extensible.size(); i++ ) {
					std::string tmp = extensible[ i ];
					if ( cur_extension_obj.find( tmp ) == cur_extension_obj.end() ) {
						skipped_fields++;
						continue;
					}
					for ( int j = 0; j < skipped_fields; j++ ) encoded += ",\n  ";
					skipped_fields = 0;
					encoded += ",\n  ";
					if ( cur_extension_obj[ tmp ].is_string() ) encoded += cur_extension_obj[ tmp ];
					else encoded += std::to_string( cur_extension_obj[ tmp ].get < double >() );
				}
			}
			encoded += ";\n\n";
		}
	}
	return encoded;
}

json IdfParser::parse_idf( std::string const & idf, size_t & index, bool & success, json const & schema ) {
	json root;
	Token token;

	while ( true ) {
		token = look_ahead( idf, index );
		if ( token == Token::END ) {
			break;
		} else if ( token == Token::NONE ) {
			success = false;
			return root;
		} else if ( token == Token::EXCLAMATION ) {
			eat_comment( idf, index );
		} else {
			std::string obj_name = parse_string( idf, index, success );
			for ( char & c : obj_name ) c = ( char ) toupper( c );
			auto tmp_umit = case_insensitive_keys.find( obj_name );
			if ( tmp_umit != case_insensitive_keys.end() ) {
				obj_name = tmp_umit->second;
			} else {
				print_out_line_error( idf, false );
				while ( token != Token::SEMICOLON && token != Token::END ) token = next_token( idf, index );
				continue;
			}
			json obj_loc = schema[ "properties" ][ obj_name ];
			json loc = obj_loc[ "legacy_idd" ];
			json obj = parse_object( idf, index, success, loc, obj_loc );
			if ( !success ) print_out_line_error( idf, true );
			std::string name = "";
			if ( !obj.is_null() ) {
				if ( obj.find( "name" ) != obj.end() ) {
					name = obj[ "name" ].get < std::string >();
					obj.erase( "name" );
				}
			}
			if ( root[ obj_name ].find( name ) != root[ obj_name ].end() ) {
				if ( name.empty() ) {
					name = obj_name + " " + std::to_string( root[ obj_name ].size() + 1 );
				} else {
					name = name + "_special_case_" + std::to_string( root[ obj_name ].size() + 1 );
				}
			}
			root[ obj_name ][ name ] = obj;
		}
	}

	return root;
}

json IdfParser::parse_object( std::string const & idf, size_t & index, bool & success,
                              json const & loc, json const & obj_loc ) {
	json root, extensible = json::object();
	json array_of_extensions = json::array();
	Token token;
	index += 1;
	int legacy_idd_index = 0, extensible_index = 0;
	success = true;
	bool was_value_parsed = false;

	while ( true ) {
		token = look_ahead( idf, index );
		if ( token == Token::NONE ) {
			success = false;
			return root;
		} else if ( token == Token::END ) {
			return root;
		} else if ( token == Token::COMMA || token == Token::SEMICOLON ) {
			if ( !was_value_parsed ) {
				int ext_size = 0;
				std::string field_name;
				if ( legacy_idd_index < loc[ "fields" ].size() ) {
					field_name = loc[ "fields" ][ legacy_idd_index ];
				} else {
					ext_size = static_cast<int>(loc[ "extensibles" ].size());
					field_name = loc[ "extensibles" ][ extensible_index % ext_size ];
					extensible_index++;
				}
				add_missing_field_value( field_name, root, extensible, obj_loc, loc, legacy_idd_index );
				if ( ext_size && extensible_index % ext_size == 0 ) {
					array_of_extensions.push_back( extensible );
					extensible.clear();
				}
			}

			legacy_idd_index++;
			was_value_parsed = false;
			next_token( idf, index );
			if ( token == Token::SEMICOLON ) {
				if ( extensible.size() ) {
					array_of_extensions.push_back( extensible );
					extensible.clear();
				}
				break;
			}
		} else if ( token == Token::EXCLAMATION ) {
			eat_comment( idf, index );
		} else if ( legacy_idd_index >= loc[ "fields" ].size() ) {
			if ( loc.find( "extensibles" ) == loc.end() ) {
				success = false;
				return root;
			}
			auto const size = loc[ "extensibles" ].size();
			std::string const field_name = loc[ "extensibles" ][ extensible_index % size ];
			auto val = parse_value( idf, index, success,
			                        obj_loc[ "patternProperties" ][ ".*" ][ "properties" ][ "extensions" ][ "items" ][ "properties" ][ field_name ] );
			extensible[ field_name ] = val;
			was_value_parsed = true;
			extensible_index++;
			if ( extensible_index && extensible_index % size == 0 ) {
				array_of_extensions.push_back( extensible );
				extensible.clear();
			}
		} else {
			was_value_parsed = true;
			std::string field = loc[ "fields" ][ legacy_idd_index ];
			auto it = obj_loc.find( "patternProperties" );
			if ( it == obj_loc.end() ) {
				if ( obj_loc.find( "properties" ) != obj_loc.end() ) {
					root[ field ] = parse_value( idf, index, success, obj_loc[ "properties" ][ field ] );
				} else {
					std::cout << "Field " << field << " was not found at line " << cur_line_num << std::endl;
				}
				legacy_idd_index++;
				continue;
			}
			auto const & tmp = obj_loc[ "patternProperties" ][ ".*" ][ "properties" ];
			if ( tmp.find( field ) == tmp.end() ) {
				if ( field == "name" ) root[ field ] = parse_string( idf, index, success );
				else std::cout << "Field " << field << " was not found at line " << cur_line_num << std::endl;
			} else {
				root[ field ] = parse_value( idf, index, success,
				                             obj_loc[ "patternProperties" ][ ".*" ][ "properties" ][ field ] );
			}
			if ( !success ) return root;
		}
	}
	if ( array_of_extensions.size() ) {
		root[ "extensions" ] = array_of_extensions;
	}
	return root;
}

void IdfParser::add_missing_field_value( std::string & field_name, json & root, json & extensible, json const & obj_loc,
                                         json const & loc, int legacy_idd_index ) {
	json tmp;
	int ext_size = 0;
	if ( obj_loc.find( "patternProperties" ) != obj_loc.end() ) {
		tmp = obj_loc[ "patternProperties" ][ ".*" ][ "properties" ];
	} else if ( obj_loc.find( "properties" ) != obj_loc.end() ) {
		tmp = obj_loc[ "properties" ][ field_name ];
	}
	if ( legacy_idd_index >= loc[ "fields" ].size() ) {
		tmp = tmp[ "extensions" ][ "items" ][ "properties" ];
		ext_size = static_cast<int>(loc[ "extensibles" ].size());
	}
	if ( tmp.find( field_name ) != tmp.end() ) {
		auto const obj_field = tmp[ field_name ];
		if ( !ext_size ) root[ field_name ] = "";
		else extensible[ field_name ] = "";
	}
}

json IdfParser::parse_number( std::string const & idf, size_t & index, bool & success ) {
	size_t save_i = index;
	eat_whitespace( idf, save_i );
	json val;
	bool is_double = false, is_sign = false, is_scientific = false;
	std::string num_str, numeric = "-+.eE0123456789";

	while ( numeric.find_first_of( idf[ save_i ] ) != std::string::npos ) {
		num_str += idf[ save_i ];
		if ( idf[ save_i ] == '.' ) {
			if ( is_double ) {
				success = false;
				return val;
			}
			is_double = true;
		} else if ( idf[ save_i ] == '-' || idf[ save_i ] == '+' ) {
			if ( is_sign && !is_scientific ) {
				success = false;
				return val;
			}
			is_sign = true;
		} else if ( idf[ save_i ] == 'e' || idf[ save_i ] == 'E' ) {
			if ( is_scientific ) {
				success = false;
				return val;
			}
			is_scientific = true;
		}
		save_i++;
	}

	assert( !num_str.empty() );
	Token token = look_ahead( idf, save_i );
	if ( token != Token::SEMICOLON && token != Token::COMMA ) {
		success = false;
		return val;
	}
	if ( is_double ) {
		val = stod( num_str, 0 );
	} else {
		try {
			val = stoi( num_str, 0 );
		} catch ( std::exception e ) {
			val = stoll( num_str, 0 );
		}
	}
	index = save_i;
	return val;
}

json IdfParser::parse_value( std::string const & idf, size_t & index, bool & success, json const & field_loc ) {
	json value;
	switch ( look_ahead( idf, index ) ) {
		case Token::STRING: {
			value = parse_string( idf, index, success );
			if ( field_loc.find( "enum" ) != field_loc.end() ) {
				for ( auto & s : field_loc[ "enum" ] ) {
					if ( icompare( s, value.get < std::string >() ) ) {
						value = s;
						break;
					}
				}
			} else if ( icompare( value.get < std::string >(), "Autosize" ) ||
			            icompare( value.get < std::string >(), "Autocalculate" ) ) {
				value = field_loc[ "anyOf" ][ 1 ][ "enum" ][ 0 ];
			}
			return value;
		}
		case Token::NUMBER: {
			size_t save_line_index = index_into_cur_line;
			size_t save_line_num = cur_line_num;
			value = parse_number( idf, index, success );
			if ( !success ) {
				cur_line_num = save_line_num;
				index_into_cur_line = save_line_index;
				success = true;
				value = parse_string( idf, index, success );
				return value;
			}
			return value;
		}
		case Token::NONE:
		case Token::END:
		case Token::EXCLAMATION:
		case Token::COMMA:
		case Token::SEMICOLON:
			break;
	}
	success = false;
	return value;
}

std::string IdfParser::parse_string( std::string const & idf, size_t & index, bool & success ) {
	eat_whitespace( idf, index );

	std::string s;
	char c;
	bool complete = false;

	while ( !complete ) {
		if ( index == idf.size() ) {
			complete = true;
			break;
		}

		c = idf[ index ];
		increment_both_index( index, index_into_cur_line );
		if ( c == ',' ) {
			complete = true;
			decrement_both_index( index, index_into_cur_line );
			break;
		} else if ( c == ';' ) {
			complete = true;
			decrement_both_index( index, index_into_cur_line );
			break;
		} else if ( c == '!' ) {
			complete = true;
			decrement_both_index( index, index_into_cur_line );
			break;
		} else if ( c == '\\' ) {
			if ( index == idf.size() ) break;
			char next_c = idf[ index ];
			increment_both_index( index, index_into_cur_line );
			if ( next_c == '"' ) {
				s += '"';
			} else if ( next_c == '\\' ) {
				s += '\\';
			} else if ( next_c == '/' ) {
				s += '/';
			} else if ( next_c == 'b' ) {
				s += '\b';
			} else if ( next_c == 't' ) {
				s += '\t';
			} else if ( next_c == 'n' ) {
				complete = false;
				break;
			} else if ( next_c == 'r' ) {
				complete = false;
				break;
			} else {
				s += c;
				s += next_c;
			}
		} else {
			s += c;
		}
	}

	if ( !complete ) {
		success = false;
		return std::string();
	}
	return rtrim( s );
}

void IdfParser::increment_both_index( size_t & index, size_t & line_index ) {
	index++;
	line_index++;
}

void IdfParser::decrement_both_index( size_t & index, size_t & line_index ) {
	index--;
	line_index--;
}

void IdfParser::print_out_line_error( std::string const & idf, bool obj_found ) {
	if ( obj_found ) std::cout << "error: \"extra field(s)\" ";
	else std::cout << "error: \"obj not found in schema\" ";
	std::cout << "at line number " << cur_line_num << " (index " << index_into_cur_line << ")\n";
	std::cout << "Line: ";
	while ( idf[ beginning_of_line_index++ ] != '\n' ) std::cout << idf[ beginning_of_line_index ];
	std::cout << std::endl;
}

void IdfParser::eat_whitespace( std::string const & idf, size_t & index ) {
	while ( index < idf.size() ) {
		switch ( idf[ index ] ) {
			case ' ':
			case '\r':
			case '\t':
				increment_both_index( index, index_into_cur_line );
				continue;
			case '\n':
				increment_both_index( index, cur_line_num );
				beginning_of_line_index = index;
				index_into_cur_line = 0;
				continue;
			default:
				return;
		}
	}
}

void IdfParser::eat_comment( std::string const & idf, size_t & index ) {
	while ( true ) {
		if ( index == idf.size() ) break;
		if ( idf[ index ] == '\n' ) {
			increment_both_index( index, cur_line_num );
			index_into_cur_line = 0;
			beginning_of_line_index = index;
			break;
		}
		increment_both_index( index, index_into_cur_line );
	}
}

IdfParser::Token IdfParser::look_ahead( std::string const & idf, size_t index ) {
	size_t save_index = index;
	size_t save_line_num = cur_line_num;
	size_t save_line_index = index_into_cur_line;
	Token token = next_token( idf, save_index );
	cur_line_num = save_line_num;
	index_into_cur_line = save_line_index;
	return token;
}

IdfParser::Token IdfParser::next_token( std::string const & idf, size_t & index ) {
	eat_whitespace( idf, index );

	if ( index == idf.size() ) {
		return Token::END;
	}

	char const c = idf[ index ];
	increment_both_index( index, index_into_cur_line );
	switch ( c ) {
		case '!':
			return Token::EXCLAMATION;
		case ',':
			return Token::COMMA;
		case ';':
			return Token::SEMICOLON;
		default:
			static std::string const search_chars( "-:.#/\\[]{}_@$%^&*()|+=<>?'\"~" );
			static std::string const numeric( ".-+0123456789" );
			if ( numeric.find_first_of( c ) != std::string::npos ) {
				return Token::NUMBER;
			} else if ( isalnum( c ) || ( std::string::npos != search_chars.find_first_of( c ) ) ) {
				return Token::STRING;
			}
			break;
	}
	decrement_both_index( index, index_into_cur_line );
	return Token::NONE;
}

void State::initialize( json & parsed_schema ) {
	stack.clear();
	schema = parsed_schema;
	stack.push_back( schema );
	json & loc = stack.back()[ "required" ];
	for ( auto & s : loc ) root_required.emplace( s.get < std::string >(), false );
}

void State::traverse( json::parse_event_t & event, json & parsed, unsigned line_num, unsigned line_index ) {
	switch ( event ) {
		case json::parse_event_t::object_start: {
			if ( is_in_extensibles or stack.back().find( "patternProperties" ) == stack.back().end() ) {
				if ( stack.back().find( "properties" ) != stack.back().end() )
					stack.push_back( stack.back()[ "properties" ] );
			} else {
				stack.push_back( stack.back()[ "patternProperties" ][ ".*" ] );
				if ( stack.back().find( "required" ) != stack.back().end() ) {
					auto & loc = stack.back()[ "required" ];
					obj_required.clear();
					for ( auto & s : loc ) obj_required.emplace( s.get < std::string >(), false );
				}
			}
			last_seen_event = event;
			break;
		}

		case json::parse_event_t::value: {
			validate( parsed, line_num, line_index );
			if ( does_key_exist ) stack.pop_back();
			does_key_exist = true;
			last_seen_event = event;
			break;
		}

		case json::parse_event_t::key: {
			std::string key = parsed;
			prev_line_index = line_index;
			prev_key_len = ( unsigned ) key.size() + 3;
			if ( need_new_object_name ) {
				cur_obj_name = key;
				cur_obj_count = 0;
				need_new_object_name = false;
				if ( cur_obj_name.find( "Parametric:" ) != std::string::npos ) {
					errors.push_back( "You must run Parametric Preprocesor for \"" + cur_obj_name + "\" at line " +
					                  std::to_string( line_num + 1 ) );
				} else if ( cur_obj_name.find( "Template" ) != std::string::npos ) {
					errors.push_back( "You must run the ExpandObjects program for \"" + cur_obj_name + "\" at line " +
					                  std::to_string( line_num + 1 ) );
				}
			}

			if ( stack.back().find( "properties" ) == stack.back().end() and key != "" ) {
				if ( stack.back().find( key ) != stack.back().end() ) {
					stack.push_back( stack.back()[ key ] );
				} else {
					errors.push_back( "Key \"" + key + "\" in object \"" + cur_obj_name + "\" at line "
					                  + std::to_string( line_num ) + " (index " + std::to_string( line_index ) +
					                  ") not found in schema" );
					does_key_exist = false;
				}
			}


			if ( !is_in_extensibles ) {
				auto req = obj_required.find( key );
				if ( req != obj_required.end() )
					req->second = true; // required field is now accounted for, for this specific object
				req = root_required.find( key );
				if ( req != root_required.end() ) req->second = true; // root_required field is now accounted for
			} else {
				auto req = extensible_required.find( key );
				if ( req != extensible_required.end() ) req->second = true;
			}

			last_seen_event = event;
			break;
		}

		case json::parse_event_t::array_start: {
			stack.push_back( stack.back()[ "items" ] );
			if ( stack.back().find( "required" ) != stack.back().end() ) {
				auto & loc = stack.back()[ "required" ];
				extensible_required.clear();
				for ( auto & s : loc ) extensible_required.emplace( s.get < std::string >(), false );
			}
			is_in_extensibles = true;
			last_seen_event = event;
			break;
		}

		case json::parse_event_t::array_end: {
			stack.pop_back();
			stack.pop_back();
			is_in_extensibles = false;
			last_seen_event = event;
			break;
		}

		case json::parse_event_t::object_end: {
			if ( is_in_extensibles ) {
				for ( auto & it : extensible_required ) {
					if ( !it.second ) {
						errors.push_back(
						"Required extensible field \"" + it.first + "\" in object \"" + cur_obj_name
						+ "\" ending at line " + std::to_string( line_num ) + " (index "
						+ std::to_string( line_index ) + ") was not provided" );
					}
					it.second = false;
				}
			} else if ( last_seen_event != json::parse_event_t::object_end ) {
				cur_obj_count++;
				for ( auto & it : obj_required ) {
					if ( !it.second ) {
						errors.push_back(
						"Required field \"" + it.first + "\" in object \"" + cur_obj_name
						+ "\" ending at line " + std::to_string( line_num ) + " (index "
						+ std::to_string( line_index ) + ") was not provided" );
					}
					it.second = false;
				}
			} else { // must be at the very end of an object now
				if ( cur_obj_name != "Version" ) stack.pop_back();
				const auto & loc = stack.back();
				if ( loc.find( "minProperties" ) != loc.end() &&
				     cur_obj_count < loc[ "minProperties" ].get < unsigned >() ) {
					errors.push_back(
					"minProperties for object \"" + cur_obj_name + "\" at line " + std::to_string( line_num ) +
					" was not met" );
				}
				if ( loc.find( "maxProperties" ) != loc.end() &&
				     cur_obj_count > loc[ "maxProperties" ].get < unsigned >() ) {
					errors.push_back(
					"maxProperties for object \"" + cur_obj_name + "\" at line " + std::to_string( line_num ) +
					" was exceeded" );
				}
				obj_required.clear();
				extensible_required.clear();
				need_new_object_name = true;
			}
			stack.pop_back();
			last_seen_event = event;
			break;
		}
	}
	if ( !stack.size() ) {
		for ( auto & it: root_required ) {
			if ( !it.second ) {
				errors.push_back( "Required object \"" + it.first + "\" was not provided in input file" );
			}
		}
	}
}

void State::validate( json & parsed, unsigned line_num, unsigned line_index ) {
	auto & loc = stack.back();

	if ( loc.find( "enum" ) != loc.end() ) {
		int i;
		auto const & enum_array = loc[ "enum" ];
		if ( parsed.is_string() ) {
			for ( i = 0; i < enum_array.size(); i++ ) {
				if ( icompare( enum_array[ i ], parsed.get < std::string >() ) ) {
					break;
				}
			}
			if ( i == enum_array.size() ) {
				errors.push_back( "In object \"" + cur_obj_name + "\" at line " + std::to_string( line_num )
				                  + ": \"" + parsed.get < std::string >() + "\" was not found in the enum" );
			}
		} else {
			for ( i = 0; i < enum_array.size(); i++ ) {
				if ( enum_array[ i ].get < int >() == parsed.get < int >() ) break;
			}
			if ( i == enum_array.size() ) {
				errors.push_back( "In object \"" + cur_obj_name + "\" at line " + std::to_string( line_num )
				                  + ": \"" + std::to_string( parsed.get < int >() ) + "\" was not found in the enum" );
			}
		}
	}
	else if ( parsed.is_number() ) {
		double val = parsed;
		if ( loc.find( "anyOf" ) != loc.end() ) {
			loc = loc[ "anyOf" ][ 0 ];
		}
		if ( loc.find( "minimum" ) != loc.end() ) {
			if ( loc.find( "exclusiveMinimum" ) != loc.end() && val <= loc[ "minimum" ].get < double >() ) {
				add_error( "exmin", val, line_num, prev_line_index + prev_key_len );
			} else if ( val < loc[ "minimum" ].get < double >() ) {
				add_error( "min", val, line_num, prev_line_index + prev_key_len );
			}
		}
		if ( loc.find( "maximum" ) != loc.end() ) {
			if ( loc.find( "exclusiveMaximum" ) != loc.end() && val >= loc[ "maximum" ].get < double >() ) {
				add_error( "exmax", val, line_num, prev_line_index + prev_key_len );
			} else if ( val > loc[ "maximum" ].get < double >() ) {
				add_error( "max", val, line_num, prev_line_index + prev_key_len );
			}
		}
		if ( loc.find( "type" ) != loc.end() && loc[ "type" ] != "number" ) {
			warnings.push_back( "In object \"" + cur_obj_name + "\" at line " + std::to_string( line_num )
			                    + ", type == " + loc[ "type" ].get < std::string >()
			                    + " but parsed value = " + std::to_string( val ) );
		}
	}
	else if ( parsed.is_string() ) {
		if ( loc.find( "anyOf" ) != loc.end() ) {
			int i;
			for ( i = 0; i < loc[ "anyOf" ].size(); i++ ) {
				if ( loc[ "anyOf" ][ i ].find( "type" ) != loc[ "anyOf" ][ i ].end() &&
				     loc[ "anyOf" ][ i ][ "type" ] == "string" )
					break;
			}
			if ( i == loc[ "anyOf" ].size() ) {
				warnings.push_back( "type == string was not found in anyOf in object \"" + cur_obj_name
				                    + "\" at line " + std::to_string( line_num ) );
			}
		}
		else {
			if ( loc.find( "type" ) != loc.end() && loc[ "type" ] != "string" ) {
				errors.push_back( "In object \"" + cur_obj_name + "\", at line " + std::to_string( line_num ) +
				                  ": type needs to be string" );
			}
		}
	}
}


namespace EnergyPlus {
// Module containing the input processor routines

// MODULE INFORMATION:
//       AUTHOR         Linda K. Lawrie
//       DATE WRITTEN   August 1997
//       MODIFIED       na
//       RE-ENGINEERED  na

// PURPOSE OF THIS MODULE:
// To provide the capabilities of reading the input data dictionary
// and input file and supplying the simulation routines with the data
// contained therein.

// METHODOLOGY EMPLOYED:

// REFERENCES:
// The input syntax is designed to allow for future flexibility without
// necessitating massive (or any) changes to this code.  Two files are
// used as key elements: (1) the input data dictionary will specify the
// sections and objects that will be allowed in the actual simulation
// input file and (2) the simulation input data file will be processed
// with the data therein being supplied to the actual simulation routines.

// OTHER NOTES:

// USE STATEMENTS:
// Use statements for data only modules
// Using/Aliasing
	using namespace DataPrecisionGlobals;
	using namespace DataStringGlobals;
	using DataGlobals::MaxNameLength;
	using DataGlobals::AutoCalculate;
	using DataGlobals::rTinyValue;
	using DataGlobals::DisplayAllWarnings;
	using DataGlobals::DisplayUnusedObjects;
	using DataGlobals::CacheIPErrorFile;
	using DataGlobals::DoingInputProcessing;
	using DataSizing::AutoSize;
	using namespace DataIPShortCuts;
	using DataSystemVariables::SortedIDD;
	using DataSystemVariables::iASCII_CR;
	using DataSystemVariables::iUnicode_end;
	using DataGlobals::DisplayInputInAudit;

	static std::string const BlankString;
	static gio::Fmt fmtLD( "*" );
	static gio::Fmt fmtA( "(A)" );
	int EchoInputFile( 0 ); // Unit number of the file echoing the IDD and input records (eplusout.audit)

// MODULE SUBROUTINES:

// Functions

// Clears the global data in InputProcessor.
// Needed for unit tests, should not be normally called.

	void
	InputProcessor::clear_state() {
		state.errors.clear();
		state.warnings.clear();
		jdf.clear();
		EchoInputFile = 0;
//	echo_stream = nullptr;
	}

	void
	InputProcessor::ProcessInput() {
		int write_stat;
		int read_stat;
		bool FileExists;

		EchoInputFile = GetNewUnitNumber();
		{
			IOFlags flags;
			flags.ACTION( "write" );
			gio::open( EchoInputFile, outputAuditFileName, flags );
			write_stat = flags.ios();
		}
		if ( write_stat != 0 ) {
			DisplayString( "Could not open (write) " + outputAuditFileName + " ." );
			ShowFatalError( "ProcessInput: Could not open file " + outputAuditFileName + " for output (write)." );
		}
		echo_stream = gio::out_stream( EchoInputFile );

		{
			IOFlags flags;
			gio::inquire( outputIperrFileName, flags );
			FileExists = flags.exists();
		}
		if ( FileExists ) {
			CacheIPErrorFile = GetNewUnitNumber();
			{
				IOFlags flags;
				flags.ACTION( "read" );
				gio::open( CacheIPErrorFile, outputIperrFileName, flags );
				read_stat = flags.ios();
			}
			if ( read_stat != 0 ) {
				ShowFatalError( "EnergyPlus: Could not open file " + outputIperrFileName + " for input (read)." );
			}
			{
				IOFlags flags;
				flags.DISPOSE( "delete" );
				gio::close( CacheIPErrorFile, flags );
			}
		}
		CacheIPErrorFile = GetNewUnitNumber();
		{
			IOFlags flags;
			flags.ACTION( "write" );
			gio::open( CacheIPErrorFile, outputIperrFileName, flags );
			write_stat = flags.ios();
		}
		if ( write_stat != 0 ) {
			DisplayString( "Could not open (write) " + outputIperrFileName );
			ShowFatalError( "ProcessInput: Could not open file " + outputIperrFileName + " for output (write)." );
		}
	}

	int
	EnergyPlus::InputProcessor::GetNumSectionsFound( std::string const & SectionWord ) {
		// PURPOSE OF THIS SUBROUTINE:
		// This function returns the number of a particular section (in input data file)
		// found in the current run.  If it can't find the section in list
		// of sections, a -1 will be returned.

		// METHODOLOGY EMPLOYED:
		// Look up section in list of sections.  If there, return the
		// number of sections of that kind found in the current input.  If not, return -1.
		if ( jdf.find( "SectionWord" ) == jdf.end() ) return -1;
		int num_sections_found = 0;
		json obj = jdf[ "SectionWord" ];
		for ( auto it = obj.begin(); it != obj.end(); ++it ) num_sections_found++;
		return num_sections_found;
	}

	int
	EnergyPlus::InputProcessor::GetNumObjectsFound( std::string const & ObjectWord ) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This function returns the number of objects (in input data file)
		// found in the current run.  If it can't find the object in list
		// of objects, a 0 will be returned.

		// METHODOLOGY EMPLOYED:
		// Look up object in list of objects.  If there, return the
		// number of objects found in the current input.  If not, return 0.

		if ( jdf.find( ObjectWord ) == jdf.end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( ObjectWord ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end()
			     || jdf.find( tmp_umit->second ) == jdf.end() ) {
				return 0;
			}
			return static_cast<int>(jdf[ tmp_umit->second ].size());
		} else {
			return static_cast<int>(jdf[ ObjectWord ].size());
		}

		if ( schema[ "properties" ].find( ObjectWord ) == schema[ "properties" ].end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( ObjectWord ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end() ) {
				ShowWarningError( "Requested Object not found in Definitions: " + ObjectWord );
			}
		}
		return 0;
	}

	void
	EnergyPlus::InputProcessor::GetObjectItem(
	std::string const & Object,
	int const Number,
	Array1S_string Alphas,
	int & NumAlphas,
	Array1S < Real64 > Numbers,
	int & NumNumbers,
	int & Status,
	Optional < Array1_bool > NumBlank,
	Optional < Array1_bool > AlphaBlank,
	Optional < Array1_string > AlphaFieldNames,
	Optional < Array1_string > NumericFieldNames
	) {

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine gets the 'number' 'object' from the IDFRecord data structure.

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int Count;
		int LoopIndex;
		std::string ObjectWord;
		std::string UCObject;
		int MaxAlphas;
		int MaxNumbers;
		int Found;
		int StartRecord;
		std::string cfld1, cfld2;
		bool GoodItem;

		json * object_in_jdf;
		if ( jdf.find( Object ) == jdf.end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( Object ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end()
			     || jdf.find( tmp_umit->second ) == jdf.end() ) {
				return;
			}
			object_in_jdf = &jdf[ tmp_umit->second ];
		} else {
			object_in_jdf = &jdf[ Object ];
		}

		json * object_in_schema;
		if ( schema[ "properties" ].find( Object ) == schema[ "properties" ].end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( Object ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end() ) {
				return;
			}
			object_in_schema = &schema[ "properties" ][ tmp_umit->second ];
		} else {
			object_in_schema = &schema[ "properties" ][ Object ];
		}

		//Autodesk:Uninit Initialize variables used uninitialized
		NumAlphas = 0; //Autodesk:Uninit Force default initialization
		NumNumbers = 0; //Autodesk:Uninit Force default initialization
		MaxAlphas = isize( Alphas, 1 );
		MaxNumbers = isize( Numbers, 1 );
		GoodItem = false;

		Status = -1;
		UCObject = MakeUPPERCase( Object );
		Alphas = "";
		Numbers = 0;

		auto obj = object_in_jdf->begin() + Number - 1;
		auto const obj_val = obj.value();
		auto const & alphas_fields = object_in_schema->at( "legacy_idd" )[ "alphas" ][ "fields" ];
		for ( int i = 0; i < alphas_fields.size(); ++i ) {
			std::string const field = alphas_fields[ i ];
			if ( field == "name" ) {
				std::string name = obj.key();
				if ( obj.key().find( "_special_case_" ) != std::string::npos ) {
					name = obj.key().substr( 0, obj.key().size() -
					                            ( obj.key().size() - obj.key().find_first_of( "_special_case_" ) ) );
				}
				Alphas( i + 1 ) = MakeUPPERCase( name );
				if ( present( AlphaBlank ) )
					AlphaBlank()( i +
					              1 ) = obj.key().empty();  // TODO is this what should be done? unit tests for the original input processor are needed
				if ( present( AlphaFieldNames ) ) AlphaFieldNames()( i + 1 ) = field;
				NumAlphas++;
				continue;
			}
			auto it = obj.value().find( field );
			if ( it != obj.value().end() ) {
				std::string val;
				if ( it.value().is_string() ) {
					json schema_obj;
					if ( object_in_schema->find( "patternProperties" ) != object_in_schema->end() ) {
						schema_obj = object_in_schema->at( "patternProperties" )[ ".*" ][ "properties" ][ field ];
					} else {
						schema_obj = object_in_schema->at( "properties" )[ field ];
					}
					if ( it.value().get < std::string >().empty() &&
					     schema_obj.find( "default" ) != schema_obj.end() ) {
						auto const & default_val = schema_obj[ "default" ];
						if ( default_val.is_string() ) {
							val = default_val.get < std::string >();
						} else {
							val = std::to_string( default_val.get < double >() );
						}
						if ( present( AlphaBlank ) ) AlphaBlank()( i + 1 ) = true;
					} else {
						val = it.value().get < std::string >();
						if ( present( AlphaBlank ) ) AlphaBlank()( i + 1 ) = val.empty();
					}
					Alphas( i + 1 ) = MakeUPPERCase( val );
				} else {
					std::string num_val = std::to_string( it.value().get < double >() );
					bool is_double = false;
					for ( auto j = num_val.find_first_of( '.' ) + 1; j != num_val.size(); j++ ) {
						if ( num_val[ j ] != '0' ) {
							is_double = true;
							break;
						}
					}
					if ( !is_double ) {
						num_val = std::to_string( it.value().get < int >() );
					}
					Alphas( i + 1 ) = num_val;
					if ( present( AlphaBlank ) ) AlphaBlank()( i + 1 ) = false;
				}
				NumAlphas++;
			} else {
				Alphas( i + 1 ) = "";
				if ( present( AlphaBlank ) ) AlphaBlank()( i + 1 ) = true;
			}
			if ( present( AlphaFieldNames ) ) AlphaFieldNames()( i + 1 ) = field;
		}

		if ( object_in_schema->at( "legacy_idd" )[ "alphas" ].find( "extensions" ) !=
		     object_in_schema->at( "legacy_idd" )[ "alphas" ].end() ) {
			auto const & alphas_extensions = object_in_schema->at( "legacy_idd" )[ "alphas" ][ "extensions" ];
			auto const extensions = obj.value()[ "extensions" ];
			int alphas_index = alphas_fields.size();
			for ( auto it = extensions.begin(); it != extensions.end(); ++it ) {
				auto const extension_obj = it.value();
				for ( auto i = 0; i < alphas_extensions.size(); i++ ) {
					std::string const field = alphas_extensions[ i ];
					if ( extension_obj.find( field ) != extension_obj.end() ) {
						if ( extension_obj[ field ].is_string() ) {
							auto const & schema_obj = object_in_schema->at(
							"patternProperties" )[ ".*" ][ "extensions" ][ field ];
							std::string val = extension_obj[ field ];
							if ( val.empty() and schema_obj.find( "default" ) != schema_obj.end() ) {
								auto const & default_val = schema_obj[ "default" ];
								if ( default_val.is_string() ) {
									val = default_val.get < std::string >();
								} else {
									val = std::to_string( default_val.get < double >() );
								}
								if ( present( AlphaBlank ) ) AlphaBlank()( alphas_index + 1 ) = true;
							} else {
								if ( present( AlphaBlank ) ) AlphaBlank()( alphas_index + 1 ) = val.empty();
							}
							Alphas( alphas_index + 1 ) = MakeUPPERCase( val );
						} else {
							double val = extension_obj[ field ];
							Alphas( alphas_index + 1 ) = std::to_string( val );
							if ( present( AlphaBlank ) ) AlphaBlank()( alphas_index + 1 ) = false;
						}
						NumAlphas++;
					} else {
						Alphas( alphas_index + 1 ) = "";
						if ( present( AlphaBlank ) ) AlphaBlank()( alphas_index + 1 ) = true;
					}
					if ( present( AlphaFieldNames ) ) AlphaFieldNames()( alphas_index + 1 ) = field;
					alphas_index++;
				}
			}
		}

		auto const & numerics_fields = object_in_schema->at( "legacy_idd" )[ "numerics" ][ "fields" ];
		for ( int i = 0; i < numerics_fields.size(); ++i ) {
			std::string const field = numerics_fields[ i ];
			auto it = obj.value().find( field );
			if ( it != obj.value().end() ) {
				if ( !it.value().is_string() ) {
					Numbers( i + 1 ) = it.value().get < double >();
					if ( present( NumBlank ) ) NumBlank()( i + 1 ) = false;
				} else {
					if ( it.value().get < std::string >().empty() ) {
						auto const & schema_obj = object_in_schema->at(
						"patternProperties" )[ ".*" ][ "properties" ][ field ]; // TODO make this account for special casing of Version like objects
						if ( schema_obj.find( "default" ) != schema_obj.end() ) {
							if ( schema_obj[ "default" ].is_string() ) Numbers( i + 1 ) = -99999;
							else Numbers( i + 1 ) = schema_obj[ "default" ].get < double >();
						} else {
							Numbers( i + 1 ) = 0;
						}
					} else {
						Numbers( i + 1 ) = -99999; // autosize and autocalculate
					}
					if ( present( NumBlank ) ) NumBlank()( i + 1 ) = it.value().get < std::string >().empty();
				}
				NumNumbers++;
			} else {
				Numbers( i + 1 ) = 0;
				if ( present( NumBlank ) )
					NumBlank()( i + 1 ) = true;
			}
			if ( present( NumericFieldNames ) ) NumericFieldNames()( i + 1 ) = field;
		}


		if ( object_in_schema->at( "legacy_idd" )[ "numerics" ].find( "extensions" ) !=
		     object_in_schema->at( "legacy_idd" )[ "numerics" ].end() ) {
			auto const & numerics_extensions = object_in_schema->at( "legacy_idd" )[ "numerics" ][ "extensions" ];
			auto const extensions = obj.value()[ "extensions" ];
			int numerics_index = numerics_fields.size();
			for ( auto it = extensions.begin(); it != extensions.end(); ++it ) {
				auto const extension_obj = it.value();
				for ( auto i = 0; i < numerics_extensions.size(); i++ ) {
					std::string const field = numerics_extensions[ i ];
					if ( extension_obj.find( field ) != extension_obj.end() ) {
						auto const val = extension_obj[ field ];
						if ( !val.is_string() ) {
							Numbers( numerics_index + 1 ) = val.get < double >();
							if ( present( NumBlank ) ) NumBlank()( numerics_index + 1 ) = false;
						} else {
							if ( val.get < std::string >().empty() ) {
								auto const & schema_obj = object_in_schema->at(
								"patternProperties" )[ ".*" ][ "properties" ][ "extensions" ][ "items" ][ "properties" ][ field ];
								if ( schema_obj.find( "default" ) != schema_obj.end() ) {
									if ( schema_obj[ "default" ].is_string() ) Numbers( numerics_index + 1 ) = -99999;
									else Numbers( numerics_index + 1 ) = schema_obj[ "default" ].get < double >();
								} else {
									Numbers( numerics_index + 1 ) = 0;
								}
							} else {
								Numbers( numerics_index + 1 ) = -99999; // default autosize value
							}
							if ( present( NumBlank ) )
								NumBlank()( numerics_index + 1 ) = val.get < std::string >().empty();
						}
						NumNumbers++;
					} else {
						auto const pattern_props = object_in_schema->find( "patternProperties" );
						if ( pattern_props != object_in_schema->end() ) {
							auto const field_in_schema = pattern_props.value()[ ".*" ][ "properties" ];
							if ( field_in_schema.find( field ) != field_in_schema.end() ) {
								if ( field_in_schema[ field ].find( "default" ) != field_in_schema[ field ].end() ) {
									auto const default_val = field_in_schema[ field ][ "default" ];
									if ( !default_val.is_string() )
										Numbers( numerics_index + 1 ) = default_val.get < double >();
									else {
										if ( it.value().get < std::string >().empty() ) {
											Numbers( numerics_index + 1 ) = 0;
										} else {
											Numbers( numerics_index + 1 ) = -99999; // autosize and autocalculate
										}
									}
								} else {
									Numbers( numerics_index + 1 ) = 0;
								}
							} else {
								std::cout << "field " << field << " not found in object " << Object << std::endl;
							}
						} else if ( object_in_schema->find( "properties" ) != object_in_schema->end() ) {
							if ( object_in_schema->at( "properties" )[ field ].find( "default" ) !=
							     object_in_schema->at( "properties" )[ field ].end() ) {
								Numbers( numerics_index + 1 ) = object_in_schema->at(
								"properties" )[ field ][ "default" ].get < double >();
							} else {
								Numbers( numerics_index + 1 ) = 0;
							}
						}
						if ( present( NumBlank ) )
							NumBlank()(
							numerics_index + 1 ) = true;
					}
					if ( present( NumericFieldNames ) ) NumericFieldNames()( numerics_index + 1 ) = field;
					numerics_index++;
				}
			}
		}
		Status = 1;
	}

	int
	EnergyPlus::InputProcessor::GetObjectItemNum(
	std::string const & ObjType, // Object Type (ref: IDD Objects)
	std::string const & ObjName // Name of the object type
	) {
		// PURPOSE OF THIS SUBROUTINE:
		// Get the occurrence number of an object of type ObjType and name ObjName

		json * obj;

		if ( jdf.find( ObjType ) == jdf.end() || jdf[ ObjType ].find( ObjName ) == jdf[ ObjType ].end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( ObjType ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end() ) {
				return -1;
			}
			obj = &jdf[ tmp_umit->second ];
		} else {
			obj = &jdf[ ObjType ];
		}

		int object_item_num = 1;
		bool found = false;
		for ( auto it = obj->begin(); it != obj->end(); ++it ) {
			if ( MakeUPPERCase( it.key() ) == MakeUPPERCase( ObjName ) ) {
				found = true;
				break;
			}
			object_item_num++;
		}

		if ( !found ) return -1;
		return object_item_num;
	}

	Real64
	EnergyPlus::InputProcessor::ProcessNumber(
	std::string const & String,
	bool & ErrorFlag
	) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function processes a string that should be numeric and
		// returns the real value of the string.

		// METHODOLOGY EMPLOYED:
		// FUNCTION ProcessNumber translates the argument (a string)
		// into a real number.  The string should consist of all
		// numeric characters (except a decimal point).  Numerics
		// with exponentiation (i.e. 1.2345E+03) are allowed but if
		// it is not a valid number an error message along with the
		// string causing the error is printed out and 0.0 is returned
		// as the value.

		// REFERENCES:
		// List directed Fortran input/output.

		static std::string const ValidNumerics(
		"0123456789.+-EeDd" ); // This had a trailing tab character: Not sure why

		Real64 rProcessNumber = 0.0;
		//  Make sure the string has all what we think numerics should have
		std::string const PString( stripped( String ) );
		std::string::size_type const StringLen( PString.length() );
		ErrorFlag = false;
		if ( StringLen == 0 ) return rProcessNumber;
		int IoStatus( 0 );
		if ( PString.find_first_not_of( ValidNumerics ) == std::string::npos ) {
			{
				IOFlags flags;
				gio::read( PString, fmtLD, flags ) >> rProcessNumber;
				IoStatus = flags.ios();
			}
			ErrorFlag = false;
		} else {
			rProcessNumber = 0.0;
			ErrorFlag = true;
		}
		if ( IoStatus != 0 ) {
			rProcessNumber = 0.0;
			ErrorFlag = true;
		}

		return rProcessNumber;

	}

	int
	EnergyPlus::InputProcessor::FindItemInList(
	std::string const & String,
	Array1_string const & ListOfItems,
	int const NumItems
	) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function looks up a string in a similar list of
		// items and returns the index of the item in the list, if
		// found.  This routine is not case insensitive and doesn't need
		// for most inputs -- they are automatically turned to UPPERCASE.
		// If you need case insensitivity use FindItem.

		for ( int Count = 1; Count <= NumItems; ++Count ) {
			if ( String == ListOfItems( Count ) ) return Count;
		}
		return 0; // Not found
	}

	int
	EnergyPlus::InputProcessor::FindItemInList(
	std::string const & String,
	Array1S_string const ListOfItems,
	int const NumItems
	) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function looks up a string in a similar list of
		// items and returns the index of the item in the list, if
		// found.  This routine is not case insensitive and doesn't need
		// for most inputs -- they are automatically turned to UPPERCASE.
		// If you need case insensitivity use FindItem.

		for ( int Count = 1; Count <= NumItems; ++Count ) {
			if ( String == ListOfItems( Count ) ) return Count;
		}
		return 0; // Not found
	}

	int
	EnergyPlus::InputProcessor::FindItemInSortedList(
	std::string const & String,
	Array1S_string const ListOfItems,
	int const NumItems
	) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function looks up a string in a similar list of
		// items and returns the index of the item in the list, if
		// found.  This routine is case insensitive.

		int Probe( 0 );
		int LBnd( 0 );
		int UBnd( NumItems + 1 );
		bool Found( false );
		while ( ( !Found ) || ( Probe != 0 ) ) {
			Probe = ( UBnd - LBnd ) / 2;
			if ( Probe == 0 ) break;
			Probe += LBnd;
			if ( equali( String, ListOfItems( Probe ) ) ) {
				Found = true;
				break;
			} else if ( lessthani( String, ListOfItems( Probe ) ) ) {
				UBnd = Probe;
			} else {
				LBnd = Probe;
			}
		}
		return Probe;
	}

	int
	EnergyPlus::InputProcessor::FindItem(
	std::string const & String,
	Array1D_string const & ListOfItems,
	int const NumItems
	) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   April 1999
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function looks up a string in a similar list of
		// items and returns the index of the item in the list, if
		// found.  This routine is case insensitive.

		int FindItem = FindItemInList( String, ListOfItems, NumItems );
		if ( FindItem != 0 ) return FindItem;

		for ( int Count = 1; Count <= NumItems; ++Count ) {
			if ( equali( String, ListOfItems( Count ) ) ) return Count;
		}
		return 0; // Not found
	}

	int
	EnergyPlus::InputProcessor::FindItem(
	std::string const & String,
	Array1S_string const ListOfItems,
	int const NumItems
	) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   April 1999
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function looks up a string in a similar list of
		// items and returns the index of the item in the list, if
		// found.  This routine is case insensitive.

		int FindItem = FindItemInList( String, ListOfItems, NumItems );
		if ( FindItem != 0 ) return FindItem;

		for ( int Count = 1; Count <= NumItems; ++Count ) {
			if ( equali( String, ListOfItems( Count ) ) ) return Count;
		}
		return 0; // Not found
	}

	std::string
	EnergyPlus::InputProcessor::MakeUPPERCase( std::string const & InputString ) {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   September 1997
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This function returns the Upper Case representation of the InputString.

		// METHODOLOGY EMPLOYED:
		// Uses the Intrinsic SCAN function to scan the lowercase representation of
		// characters (DataStringGlobals) for each character in the given string.

		// Locals
		// FUNCTION ARGUMENT DEFINITIONS:
		// MaxInputLineLength because of PowerStation Compiler
		// otherwise could say (CHARACTER(len=LEN(InputString))

		std::string ResultString( InputString );

		for ( std::string::size_type i = 0, e = len( InputString ); i < e; ++i ) {
			int const curCharVal = int( InputString[ i ] );
			if ( ( 97 <= curCharVal && curCharVal <= 122 ) ||
			     ( 224 <= curCharVal && curCharVal <= 255 ) ) { // lowercase ASCII and accented characters
				ResultString[ i ] = char( curCharVal - 32 );
			}
		}

		return ResultString;

	}

	void
	EnergyPlus::InputProcessor::VerifyName(
	std::string const & NameToVerify,
	Array1D_string const & NamesList,
	int const NumOfNames,
	bool & ErrorFound,
	bool & IsBlank,
	std::string const & StringToDisplay
	) {

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   February 2000
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine verifys that a new name can be added to the
		// list of names for this item (i.e., that there isn't one of that
		// name already and that this name is not blank).

		int Found;

		ErrorFound = false;
		if ( NumOfNames > 0 ) {
			Found = FindItem( NameToVerify, NamesList, NumOfNames );
			if ( Found != 0 ) {
				ShowSevereError( StringToDisplay + ", duplicate name=" + NameToVerify );
				ErrorFound = true;
			}
		}

		if ( NameToVerify.empty() ) {
			ShowSevereError( StringToDisplay + ", cannot be blank" );
			ErrorFound = true;
			IsBlank = true;
		} else {
			IsBlank = false;
		}

	}

	void
	EnergyPlus::InputProcessor::VerifyName(
	std::string const & NameToVerify,
	Array1S_string const NamesList,
	int const NumOfNames,
	bool & ErrorFound,
	bool & IsBlank,
	std::string const & StringToDisplay
	) {

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   February 2000
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine verifys that a new name can be added to the
		// list of names for this item (i.e., that there isn't one of that
		// name already and that this name is not blank).

		int Found;

		ErrorFound = false;
		if ( NumOfNames > 0 ) {
			Found = FindItem( NameToVerify, NamesList, NumOfNames );
			if ( Found != 0 ) {
				ShowSevereError( StringToDisplay + ", duplicate name=" + NameToVerify );
				ErrorFound = true;
			}
		}

		if ( NameToVerify.empty() ) {
			ShowSevereError( StringToDisplay + ", cannot be blank" );
			ErrorFound = true;
			IsBlank = true;
		} else {
			IsBlank = false;
		}

	}

	void
	EnergyPlus::InputProcessor::RangeCheck(
	bool & ErrorsFound, // Set to true if error detected
	std::string const & WhatFieldString, // Descriptive field for string
	std::string const & WhatObjectString, // Descriptive field for object, Zone Name, etc.
	std::string const & ErrorLevel, // 'Warning','Severe','Fatal')
	Optional_string_const LowerBoundString, // String for error message, if applicable
	Optional_bool_const LowerBoundCondition, // Condition for error condition, if applicable
	Optional_string_const UpperBoundString, // String for error message, if applicable
	Optional_bool_const UpperBoundCondition, // Condition for error condition, if applicable
	Optional_string_const ValueString, // Value with digits if to be displayed with error
	Optional_string_const WhatObjectName // ObjectName -- used for error messages
	) {

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   July 2000
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is a general purpose "range check" routine for GetInput routines.
		// Using the standard "ErrorsFound" logical, this routine can produce a reasonable
		// error message to describe the situation in addition to setting the ErrorsFound variable
		// to true.

		std::string ErrorString; // Uppercase representation of ErrorLevel
		std::string Message1;
		std::string Message2;

		bool Error = false;
		if ( present( UpperBoundCondition ) ) {
			if ( !UpperBoundCondition ) Error = true;
		}
		if ( present( LowerBoundCondition ) ) {
			if ( !LowerBoundCondition ) Error = true;
		}

		if ( Error ) {
			ConvertCaseToUpper( ErrorLevel, ErrorString );
			Message1 = WhatObjectString;
			if ( present( WhatObjectName ) ) Message1 += "=\"" + WhatObjectName + "\", out of range data";
			Message2 = "Out of range value field=" + WhatFieldString;
			if ( present( ValueString ) ) Message2 += ", Value=[" + ValueString + ']';
			Message2 += ", range={";
			if ( present( LowerBoundString ) ) Message2 += LowerBoundString;
			if ( present( LowerBoundString ) && present( UpperBoundString ) ) {
				Message2 += " and " + UpperBoundString;
			} else if ( present( UpperBoundString ) ) {
				Message2 += UpperBoundString;
			}
			Message2 += "}";

			{
				auto const errorCheck( ErrorString[ 0 ] );

				if ( ( errorCheck == 'W' ) || ( errorCheck == 'w' ) ) {
					ShowWarningError( Message1 );
					ShowContinueError( Message2 );

				} else if ( ( errorCheck == 'S' ) || ( errorCheck == 's' ) ) {
					ShowSevereError( Message1 );
					ShowContinueError( Message2 );
					ErrorsFound = true;

				} else if ( ( errorCheck == 'F' ) || ( errorCheck == 'f' ) ) {
					ShowSevereError( Message1 );
					ShowContinueError( Message2 );
					ShowFatalError( "Program terminates due to preceding condition(s)." );

				} else {
					ShowSevereError( Message1 );
					ShowContinueError( Message2 );
					ErrorsFound = true;
				}
			}
		}
	}

	int
	EnergyPlus::InputProcessor::GetNumRangeCheckErrorsFound() {

		// FUNCTION INFORMATION:
		//       AUTHOR         Linda K. Lawrie
		//       DATE WRITTEN   July 2000
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// This function returns the number of OutOfRange errors found during
		// input processing.
//	return NumOutOfRangeErrorsFound;
		// TODO: Fix this
		return 0;
	}

	void
	EnergyPlus::InputProcessor::GetObjectDefMaxArgs(
	std::string const & ObjectWord, // Object for definition
	int & NumArgs, // How many arguments (max) this Object can have
	int & NumAlpha, // How many Alpha arguments (max) this Object can have
	int & NumNumeric // How many Numeric arguments (max) this Object can have
	) {
		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine returns maximum argument limits (total, alphas, numerics) of an Object from the IDD.
		// These dimensions (not sure what one can use the total for) can be used to dynamically dimension the
		// arrays in the GetInput routines.

		NumArgs = 0;
		NumAlpha = 0;
		NumNumeric = 0;
		json * object;
		if ( schema[ "properties" ].find( ObjectWord ) == schema[ "properties" ].end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( ObjectWord ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end() ) {
				ShowSevereError(
				"GetObjectDefMaxArgs: Did not find object=\"" + ObjectWord + "\" in list of objects." );
				return;
			}
			object = &schema[ "properties" ][ tmp_umit->second ];
		} else {
			object = &schema[ "properties" ][ ObjectWord ];
		}

		const json & legacy_idd = object->at( "legacy_idd" );

		json * objects;
		if ( jdf.find( ObjectWord ) == jdf.end() ) {
			auto tmp_umit = InputProcessor::idf_parser.case_insensitive_keys.find( MakeUPPERCase( ObjectWord ) );
			if ( tmp_umit == InputProcessor::idf_parser.case_insensitive_keys.end() ) {
				ShowSevereError(
				"GetObjectDefMaxArgs: Did not find object=\"" + ObjectWord + "\" in list of objects." );
				return;
			}
			objects = &jdf[ tmp_umit->second ];
		} else {
			objects = &jdf[ ObjectWord ];
		}

		size_t max_size = 0;
		for ( auto const obj : *objects ) {
			if ( obj.find( "extensions" ) != obj.end() ) {
				auto const size = obj[ "extensions" ].size();
				if ( size > max_size ) max_size = size;
			}
		}

		if ( legacy_idd.find( "alphas" ) != legacy_idd.end() ) {
			json const alphas = legacy_idd[ "alphas" ];
			if ( alphas.find( "fields" ) != alphas.end() ) {
				NumAlpha += alphas[ "fields" ].size();
			}
			if ( alphas.find( "extensions" ) != alphas.end() ) {
				NumAlpha += alphas[ "extensions" ].size() * max_size;
			}
		}
		if ( legacy_idd.find( "numerics" ) != legacy_idd.end() ) {
			json const numerics = legacy_idd[ "numerics" ];
			if ( numerics.find( "fields" ) != numerics.end() ) {
				NumNumeric += numerics[ "fields" ].size();
			}
			if ( numerics.find( "extensions" ) != numerics.end() ) {
				NumNumeric += numerics[ "extensions" ].size() * max_size;
			}
		}
		NumArgs = NumAlpha + NumNumeric;
	}

	void
	EnergyPlus::InputProcessor::PreProcessorCheck( bool & PreP_Fatal ) // True if a preprocessor flags a fatal error
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   August 2005
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This routine checks for existance of "Preprocessor Message" object and
		// performs appropriate action.

		// METHODOLOGY EMPLOYED:
		// na

		// REFERENCES:
		// Preprocessor Message,
		//    \memo This object does not come from a user input.  This is generated by a pre-processor
		//    \memo so that various conditions can be gracefully passed on by the InputProcessor.
		//    A1,        \field preprocessor name
		//    A2,        \field error severity
		//               \note Depending on type, InputProcessor may terminate the program.
		//               \type choice
		//               \key warning
		//               \key severe
		//               \key fatal
		//    A3,        \field message line 1
		//    A4,        \field message line 2
		//    A5,        \field message line 3
		//    A6,        \field message line 4
		//    A7,        \field message line 5
		//    A8,        \field message line 6
		//    A9,        \field message line 7
		//    A10,       \field message line 8
		//    A11,       \field message line 9
		//    A12;       \field message line 10

		// Using/Aliasing
		using namespace DataIPShortCuts;

		int NumAlphas; // Used to retrieve names from IDF
		int NumNumbers; // Used to retrieve rNumericArgs from IDF
		int IOStat; // Could be used in the Get Routines, not currently checked
		int NumParams; // Total Number of Parameters in 'Output:PreprocessorMessage' Object
		int NumPrePM; // Number of Preprocessor Message objects in IDF
		int CountP;
		int CountM;
		std::string Multiples;

		cCurrentModuleObject = "Output:PreprocessorMessage";
		NumPrePM = InputProcessor::GetNumObjectsFound( cCurrentModuleObject );
		if ( NumPrePM > 0 ) {
			GetObjectDefMaxArgs( cCurrentModuleObject, NumParams, NumAlphas, NumNumbers );
			cAlphaArgs( { 1, NumAlphas } ) = BlankString;
			for ( CountP = 1; CountP <= NumPrePM; ++CountP ) {
				InputProcessor::GetObjectItem( cCurrentModuleObject, CountP, cAlphaArgs, NumAlphas, rNumericArgs,
				                               NumNumbers, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks,
				                               cAlphaFieldNames, cNumericFieldNames );
				if ( cAlphaArgs( 1 ).empty() ) cAlphaArgs( 1 ) = "Unknown";
				if ( NumAlphas > 3 ) {
					Multiples = "s";
				} else {
					Multiples = BlankString;
				}
				if ( cAlphaArgs( 2 ).empty() ) cAlphaArgs( 2 ) = "Unknown";
				{
					auto const errorType( uppercased( cAlphaArgs( 2 ) ) );
					if ( errorType == "INFORMATION" ) {
						ShowMessage( cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) +
						             "\" has the following Information message" + Multiples + ':' );
					} else if ( errorType == "WARNING" ) {
						ShowWarningError( cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) +
						                  "\" has the following Warning condition" + Multiples + ':' );
					} else if ( errorType == "SEVERE" ) {
						ShowSevereError(
						cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) +
						"\" has the following Severe condition" +
						Multiples + ':' );
					} else if ( errorType == "FATAL" ) {
						ShowSevereError(
						cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) +
						"\" has the following Fatal condition" +
						Multiples + ':' );
						PreP_Fatal = true;
					} else {
						ShowSevereError(
						cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\" has the following " +
						cAlphaArgs( 2 ) +
						" condition" + Multiples + ':' );
					}
				}
				CountM = 3;
				if ( CountM > NumAlphas ) {
					ShowContinueError( cCurrentModuleObject + " was blank.  Check " + cAlphaArgs( 1 ) +
					                   " audit trail or error file for possible reasons." );
				}
				while ( CountM <= NumAlphas ) {
					if ( len( cAlphaArgs( CountM ) ) == MaxNameLength ) {
						ShowContinueError( cAlphaArgs( CountM ) + cAlphaArgs( CountM + 1 ) );
						CountM += 2;
					} else {
						ShowContinueError( cAlphaArgs( CountM ) );
						++CountM;
					}
				}
			}
		}

	}

	void
	InputProcessor::PreScanReportingVariables() {
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   July 2010

		// PURPOSE OF THIS SUBROUTINE:
		// This routine scans the input records and determines which output variables
		// are actually being requested for the run so that the OutputProcessor will only
		// consider those variables for output.  (At this time, all metered variables are
		// allowed to pass through).

		// METHODOLOGY EMPLOYED:
		// Uses internal records and structures.
		// Looks at:
		// Output:Variable
		// Meter:Custom
		// Meter:CustomDecrement
		// Meter:CustomDifference
		// Output:Table:Monthly
		// Output:Table:TimeBins
		// Output:Table:SummaryReports
		// EnergyManagementSystem:Sensor
		// EnergyManagementSystem:OutputVariable

		// Using/Aliasing
		using namespace DataOutputs;

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const OutputVariable( "Output:Variable" );
		static std::string const MeterCustom( "Meter:Custom" );
		static std::string const MeterCustomDecrement( "Meter:CustomDecrement" );
//		static std::string const MeterCustomDifference( "METER:CUSTOMDIFFERENCE" );
		static std::string const OutputTableMonthly( "Output:Table:Monthly" );
		static std::string const OutputTableAnnual( "Output:Table:Annual" );
		static std::string const OutputTableTimeBins( "Output:Table:TimeBins" );
		static std::string const OutputTableSummaries( "Output:Table:SummaryReports" );
		static std::string const EMSSensor( "EnergyManagementSystem:Sensor" );
		static std::string const EMSOutputVariable( "EnergyManagementSystem:OutputVariable" );

		OutputVariablesForSimulation.allocate( 10000 );
		MaxConsideredOutputVariables = 10000;

		// Output Variable
		auto jdf_objects = jdf.find( OutputVariable );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				if ( !fields.at( "key_value" ).empty() ) {
					InputProcessor::AddRecordToOutputVariableStructure( fields.at( "key_value" ),
					                                                    fields.at( "variable_name" ) );
				} else {
					InputProcessor::AddRecordToOutputVariableStructure( "*", fields.at( "variable_name" ) );
				}
			}
		}

		jdf_objects = jdf.find( MeterCustom );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				for ( auto const & extensions : fields[ "extensions" ] ) {
					if ( !obj.key().empty() ) {
						InputProcessor::AddRecordToOutputVariableStructure( extensions.at( "key_name" ),
						                                                    extensions.at(
						                                                    "output_variable_or_meter_name" ) );
					} else {
						InputProcessor::AddRecordToOutputVariableStructure( "*", extensions.at(
						"output_variable_or_meter_name" ) );
					}
				}
			}
		}

		jdf_objects = jdf.find( MeterCustomDecrement );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				for ( auto const & extensions : fields[ "extensions" ] ) {
					if ( !obj.key().empty() ) {
						InputProcessor::AddRecordToOutputVariableStructure( extensions.at( "key_name" ),
						                                                    extensions.at(
						                                                    "output_variable_or_meter_name" ) );
					} else {
						InputProcessor::AddRecordToOutputVariableStructure( "*", extensions.at( "output_variable_or_meter_name" ) );
					}
				}
			}
		}

		jdf_objects = jdf.find( EMSSensor );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				if ( !fields.at( "output_variable_or_output_meter_index_key_name" ).empty() ) {
					InputProcessor::AddRecordToOutputVariableStructure(
					fields.at( "output_variable_or_output_meter_index_key_name" ),
					fields.at( "output_variable_or_output_meter_name" ) );
				} else {
					InputProcessor::AddRecordToOutputVariableStructure( "*", fields.at( "output_variable_or_output_meter_name" ) );
				}
			}
		}

		jdf_objects = jdf.find( EMSOutputVariable );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				InputProcessor::AddRecordToOutputVariableStructure( "*", obj.key() );
			}
		}

		jdf_objects = jdf.find( OutputTableTimeBins );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				if ( !obj.key().empty() ) {
					InputProcessor::AddRecordToOutputVariableStructure( obj.key(), fields.at( "key_value" ) );
				} else {
					InputProcessor::AddRecordToOutputVariableStructure( "*", fields.at( "key_value" ) );
				}
			}
		}

		jdf_objects = jdf.find( OutputTableMonthly );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				for ( auto const & extensions : fields[ "extensions" ] ) {
					InputProcessor::AddRecordToOutputVariableStructure( "*", extensions.at( "variable_or_meter_name" ) );
				}
			}
		}

		jdf_objects = jdf.find( OutputTableAnnual );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				for ( auto const & extensions : fields[ "extensions" ] ) {
					InputProcessor::AddRecordToOutputVariableStructure( "*", extensions.at(
					"variable_or_meter_or_ems_variable_or_field_name" ) );
				}
			}
		}

		jdf_objects = jdf.find( OutputTableSummaries );
		if ( jdf_objects != jdf.end() ) {
			auto const & jdf_object = jdf_objects.value();
			for ( auto obj = jdf_object.begin(); obj != jdf_object.end(); ++obj ) {
				json const & fields = obj.value();
				for ( auto const & extensions : fields[ "extensions" ] ) {
					auto const report_name = MakeUPPERCase( extensions.at( "report_name" ) );
					if ( report_name == "ALLMONTHLY" || report_name == "ALLSUMMARYANDMONTHLY" ) {
						for ( int i = 1; i <= NumMonthlyReports; ++i ) {
							InputProcessor::AddVariablesForMonthlyReport( MonthlyNamedReports( i ) );
						}
					} else {
						InputProcessor::AddVariablesForMonthlyReport( report_name );
					}
				}
			}
		}

		if ( NumConsideredOutputVariables > 0 ) {
			OutputVariablesForSimulation.redimension( NumConsideredOutputVariables );
			MaxConsideredOutputVariables = NumConsideredOutputVariables;
		}
	}

	void
	InputProcessor::AddVariablesForMonthlyReport( std::string const & reportName ) {

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   July 2010
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This routine adds specific variables to the Output Variables for Simulation
		// Structure. Note that only non-metered variables need to be added here.  Metered
		// variables are automatically included in the minimized output variable structure.

		if ( reportName == "ZONECOOLINGSUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE AIR SYSTEM SENSIBLE COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR WETBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "ZONE TOTAL INTERNAL LATENT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE TOTAL INTERNAL LATENT GAIN RATE" );

		} else if ( reportName == "ZONEHEATINGSUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE AIR SYSTEM SENSIBLE HEATING ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "ZONE AIR SYSTEM SENSIBLE HEATING RATE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );

		} else if ( reportName == "ZONEELECTRICSUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE LIGHTS ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE ELECTRIC EQUIPMENT ELECTRIC ENERGY" );

		} else if ( reportName == "SPACEGAINSMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE PEOPLE TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE LIGHTS TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE ELECTRIC EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE GAS EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE HOT WATER EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE STEAM EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE OTHER EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE INFILTRATION SENSIBLE HEAT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE INFILTRATION SENSIBLE HEAT LOSS ENERGY" );

		} else if ( reportName == "PEAKSPACEGAINSMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE PEOPLE TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE LIGHTS TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE ELECTRIC EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE GAS EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE HOT WATER EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE STEAM EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE OTHER EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE INFILTRATION SENSIBLE HEAT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE INFILTRATION SENSIBLE HEAT LOSS ENERGY" );

		} else if ( reportName == "SPACEGAINCOMPONENTSATCOOLINGPEAKMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE AIR SYSTEM SENSIBLE COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "ZONE PEOPLE TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE LIGHTS TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE ELECTRIC EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE GAS EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE HOT WATER EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE STEAM EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE OTHER EQUIPMENT TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE INFILTRATION SENSIBLE HEAT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE INFILTRATION SENSIBLE HEAT LOSS ENERGY" );

		} else if ( reportName == "SETPOINTSNOTMETWITHTEMPERATURESMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE HEATING SETPOINT NOT MET TIME" );
			AddRecordToOutputVariableStructure( "*", "ZONE MEAN AIR TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "ZONE HEATING SETPOINT NOT MET WHILE OCCUPIED TIME" );
			AddRecordToOutputVariableStructure( "*", "ZONE COOLING SETPOINT NOT MET TIME" );
			AddRecordToOutputVariableStructure( "*", "ZONE COOLING SETPOINT NOT MET WHILE OCCUPIED TIME" );

		} else if ( reportName == "COMFORTREPORTSIMPLE55MONTHLY" ) {
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE THERMAL COMFORT ASHRAE 55 SIMPLE MODEL SUMMER CLOTHES NOT COMFORTABLE TIME" );
			AddRecordToOutputVariableStructure( "*", "ZONE MEAN AIR TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE THERMAL COMFORT ASHRAE 55 SIMPLE MODEL WINTER CLOTHES NOT COMFORTABLE TIME" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE THERMAL COMFORT ASHRAE 55 SIMPLE MODEL SUMMER OR WINTER CLOTHES NOT COMFORTABLE TIME" );

		} else if ( reportName == "UNGLAZEDTRANSPIREDSOLARCOLLECTORSUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SOLAR COLLECTOR SYSTEM EFFICIENCY" );
			AddRecordToOutputVariableStructure( "*", "SOLAR COLLECTOR OUTSIDE FACE SUCTION VELOCITY" );
			AddRecordToOutputVariableStructure( "*", "SOLAR COLLECTOR SENSIBLE HEATING RATE" );

		} else if ( reportName == "OCCUPANTCOMFORTDATASUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "PEOPLE OCCUPANT COUNT" );
			AddRecordToOutputVariableStructure( "*", "PEOPLE AIR TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "PEOPLE AIR RELATIVE HUMIDITY" );
			AddRecordToOutputVariableStructure( "*", "ZONE THERMAL COMFORT FANGER MODEL PMV" );
			AddRecordToOutputVariableStructure( "*", "ZONE THERMAL COMFORT FANGER MODEL PPD" );

		} else if ( reportName == "CHILLERREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "CHILLER ELECTRIC ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "CHILLER ELECTRIC POWER" );
			AddRecordToOutputVariableStructure( "*", "CHILLER EVAPORATOR COOLING ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "CHILLER CONDENSER HEAT TRANSFER ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "CHILLER COP" );

		} else if ( reportName == "TOWERREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "COOLING TOWER FAN ELECTRIC ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "COOLING TOWER FAN ELECTRIC POWER" );
			AddRecordToOutputVariableStructure( "*", "COOLING TOWER HEAT TRANSFER RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING TOWER INLET TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "COOLING TOWER OUTLET TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "COOLING TOWER MASS FLOW RATE" );

		} else if ( reportName == "BOILERREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "BOILER HEATING ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "BOILER GAS CONSUMPTION" ); // on meter
			AddRecordToOutputVariableStructure( "*", "BOILER HEATING ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "BOILER HEATING RATE" );
			AddRecordToOutputVariableStructure( "*", "BOILER GAS CONSUMPTION RATE" );
			AddRecordToOutputVariableStructure( "*", "BOILER INLET TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "BOILER OUTLET TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "BOILER MASS FLOW RATE" );
			AddRecordToOutputVariableStructure( "*", "BOILER ANCILLARY ELECTRIC POWER" );

		} else if ( reportName == "DXREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "COOLING COIL TOTAL COOLING ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "COOLING COIL ELECTRIC ENERGY" ); // on meter
			AddRecordToOutputVariableStructure( "*", "COOLING COIL SENSIBLE COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL LATENT COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL CRANKCASE HEATER ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL RUNTIME FRACTION" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL TOTAL COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL SENSIBLE COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL LATENT COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL ELECTRIC POWER" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL CRANKCASE HEATER ELECTRIC POWER" );

		} else if ( reportName == "WINDOWREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW TRANSMITTED SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW TRANSMITTED BEAM SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW TRANSMITTED DIFFUSE SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW HEAT GAIN RATE" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW HEAT LOSS RATE" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW INSIDE FACE GLAZING CONDENSATION STATUS" );
			AddRecordToOutputVariableStructure( "*", "SURFACE SHADING DEVICE IS ON TIME FRACTION" );
			AddRecordToOutputVariableStructure( "*", "SURFACE STORM WINDOW ON OFF STATUS" );

		} else if ( reportName == "WINDOWENERGYREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW TRANSMITTED SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW TRANSMITTED BEAM SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW TRANSMITTED DIFFUSE SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW HEAT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "SURFACE WINDOW HEAT LOSS ENERGY" );

		} else if ( reportName == "WINDOWZONESUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOWS TOTAL HEAT GAIN RATE" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOWS TOTAL HEAT LOSS RATE" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOWS TOTAL TRANSMITTED SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE EXTERIOR WINDOWS TOTAL TRANSMITTED BEAM SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE EXTERIOR WINDOWS TOTAL TRANSMITTED DIFFUSE SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE INTERIOR WINDOWS TOTAL TRANSMITTED DIFFUSE SOLAR RADIATION RATE" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE INTERIOR WINDOWS TOTAL TRANSMITTED BEAM SOLAR RADIATION RATE" );

		} else if ( reportName == "WINDOWENERGYZONESUMMARYMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOWS TOTAL HEAT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOWS TOTAL HEAT LOSS ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOWS TOTAL TRANSMITTED SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE EXTERIOR WINDOWS TOTAL TRANSMITTED BEAM SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE EXTERIOR WINDOWS TOTAL TRANSMITTED DIFFUSE SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE INTERIOR WINDOWS TOTAL TRANSMITTED DIFFUSE SOLAR RADIATION ENERGY" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE INTERIOR WINDOWS TOTAL TRANSMITTED BEAM SOLAR RADIATION ENERGY" );

		} else if ( reportName == "AVERAGEOUTDOORCONDITIONSMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR WETBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DEWPOINT TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE WIND SPEED" );
			AddRecordToOutputVariableStructure( "*", "SITE SKY TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE DIFFUSE SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE DIRECT SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE RAIN STATUS" );

		} else if ( reportName == "OUTDOORCONDITIONSMAXIMUMDRYBULBMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR WETBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DEWPOINT TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE WIND SPEED" );
			AddRecordToOutputVariableStructure( "*", "SITE SKY TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE DIFFUSE SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE DIRECT SOLAR RADIATION RATE PER AREA" );

		} else if ( reportName == "OUTDOORCONDITIONSMINIMUMDRYBULBMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR WETBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DEWPOINT TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE WIND SPEED" );
			AddRecordToOutputVariableStructure( "*", "SITE SKY TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE DIFFUSE SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE DIRECT SOLAR RADIATION RATE PER AREA" );

		} else if ( reportName == "OUTDOORCONDITIONSMAXIMUMWETBULBMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR WETBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DEWPOINT TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE WIND SPEED" );
			AddRecordToOutputVariableStructure( "*", "SITE SKY TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE DIFFUSE SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE DIRECT SOLAR RADIATION RATE PER AREA" );

		} else if ( reportName == "OUTDOORCONDITIONSMAXIMUMDEWPOINTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DEWPOINT TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR DRYBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE OUTDOOR AIR WETBULB TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE WIND SPEED" );
			AddRecordToOutputVariableStructure( "*", "SITE SKY TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE DIFFUSE SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE DIRECT SOLAR RADIATION RATE PER AREA" );

		} else if ( reportName == "OUTDOORGROUNDCONDITIONSMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE GROUND TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE SURFACE GROUND TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE DEEP GROUND TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE MAINS WATER TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "SITE GROUND REFLECTED SOLAR RADIATION RATE PER AREA" );
			AddRecordToOutputVariableStructure( "*", "SITE SNOW ON GROUND STATUS" );

		} else if ( reportName == "WINDOWACREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER TOTAL COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER TOTAL COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER SENSIBLE COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER LATENT COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER TOTAL COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER SENSIBLE COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER LATENT COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "ZONE WINDOW AIR CONDITIONER ELECTRIC POWER" );

		} else if ( reportName == "WATERHEATERREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "WATER HEATER TOTAL DEMAND HEAT TRANSFER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER USE SIDE HEAT TRANSFER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER BURNER HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER GAS CONSUMPTION" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER TOTAL DEMAND HEAT TRANSFER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER LOSS DEMAND ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER HEAT LOSS ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER TANK TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER HEAT RECOVERY SUPPLY ENERGY" );
			AddRecordToOutputVariableStructure( "*", "WATER HEATER SOURCE ENERGY" );

		} else if ( reportName == "GENERATORREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "GENERATOR PRODUCED ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR DIESEL CONSUMPTION" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR GAS CONSUMPTION" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR PRODUCED ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR TOTAL HEAT RECOVERY" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR JACKET HEAT RECOVERY ENERGY" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR LUBE HEAT RECOVERY" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR EXHAUST HEAT RECOVERY ENERGY" );
			AddRecordToOutputVariableStructure( "*", "GENERATOR EXHAUST AIR TEMPERATURE" );

		} else if ( reportName == "DAYLIGHTINGREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "SITE EXTERIOR BEAM NORMAL ILLUMINANCE" );
			AddRecordToOutputVariableStructure( "*", "DAYLIGHTING LIGHTING POWER MULTIPLIER" );
			AddRecordToOutputVariableStructure( "*", "DAYLIGHTING LIGHTING POWER MULTIPLIER" );
			AddRecordToOutputVariableStructure( "*", "DAYLIGHTING REFERENCE POINT 1 ILLUMINANCE" );
			AddRecordToOutputVariableStructure( "*", "DAYLIGHTING REFERENCE POINT 1 GLARE INDEX" );
			AddRecordToOutputVariableStructure( "*",
			                                    "DAYLIGHTING REFERENCE POINT 1 GLARE INDEX SETPOINT EXCEEDED TIME" );
			AddRecordToOutputVariableStructure( "*",
			                                    "DAYLIGHTING REFERENCE POINT 1 DAYLIGHT ILLUMINANCE SETPOINT EXCEEDED TIME" );
			AddRecordToOutputVariableStructure( "*", "DAYLIGHTING REFERENCE POINT 2 ILLUMINANCE" );
			AddRecordToOutputVariableStructure( "*", "DAYLIGHTING REFERENCE POINT 2 GLARE INDEX" );
			AddRecordToOutputVariableStructure( "*",
			                                    "DAYLIGHTING REFERENCE POINT 2 GLARE INDEX SETPOINT EXCEEDED TIME" );
			AddRecordToOutputVariableStructure( "*",
			                                    "DAYLIGHTING REFERENCE POINT 2 DAYLIGHT ILLUMINANCE SETPOINT EXCEEDED TIME" );

		} else if ( reportName == "COILREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "HEATING COIL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "HEATING COIL HEATING RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL SENSIBLE COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL TOTAL COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL TOTAL COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL SENSIBLE COOLING RATE" );
			AddRecordToOutputVariableStructure( "*", "COOLING COIL WETTED AREA FRACTION" );

		} else if ( reportName == "PLANTLOOPDEMANDREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "PLANT SUPPLY SIDE COOLING DEMAND RATE" );
			AddRecordToOutputVariableStructure( "*", "PLANT SUPPLY SIDE HEATING DEMAND RATE" );

		} else if ( reportName == "FANREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "FAN ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "FAN RISE IN AIR TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "FAN ELECTRIC POWER" );

		} else if ( reportName == "PUMPREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "PUMP ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "PUMP FLUID HEAT GAIN ENERGY" );
			AddRecordToOutputVariableStructure( "*", "PUMP ELECTRIC POWER" );
			AddRecordToOutputVariableStructure( "*", "PUMP SHAFT POWER" );
			AddRecordToOutputVariableStructure( "*", "PUMP FLUID HEAT GAIN RATE" );
			AddRecordToOutputVariableStructure( "*", "PUMP OUTLET TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "PUMP MASS FLOW RATE" );

		} else if ( reportName == "CONDLOOPDEMANDREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "PLANT SUPPLY SIDE COOLING DEMAND RATE" );
			AddRecordToOutputVariableStructure( "*", "PLANT SUPPLY SIDE HEATING DEMAND RATE" );
			AddRecordToOutputVariableStructure( "*", "PLANT SUPPLY SIDE INLET TEMPERATURE" );
			AddRecordToOutputVariableStructure( "*", "PLANT SUPPLY SIDE OUTLET TEMPERATURE" );

		} else if ( reportName == "ZONETEMPERATUREOSCILLATIONREPORTMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE OSCILLATING TEMPERATURES TIME" );
			AddRecordToOutputVariableStructure( "*", "ZONE PEOPLE OCCUPANT COUNT" );

		} else if ( reportName == "AIRLOOPSYSTEMENERGYANDWATERUSEMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HOT WATER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM STEAM ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM CHILLED WATER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM GAS ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM WATER VOLUME" );

		} else if ( reportName == "AIRLOOPSYSTEMCOMPONENTLOADSMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM FAN AIR HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM COOLING COIL TOTAL COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEATING COIL TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEAT EXCHANGER TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEAT EXCHANGER TOTAL COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HUMIDIFIER TOTAL HEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM EVAPORATIVE COOLER TOTAL COOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM DESICCANT DEHUMIDIFIER TOTAL COOLING ENERGY" );

		} else if ( reportName == "AIRLOOPSYSTEMCOMPONENTENERGYUSEMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM FAN ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEATING COIL HOT WATER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM COOLING COIL CHILLED WATER ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM DX HEATING COIL ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM DX COOLING COIL ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEATING COIL ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEATING COIL GAS ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HEATING COIL STEAM ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM HUMIDIFIER ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM EVAPORATIVE COOLER ELECTRIC ENERGY" );
			AddRecordToOutputVariableStructure( "*", "AIR SYSTEM DESICCANT DEHUMIDIFIER ELECTRIC ENERGY" );

		} else if ( reportName == "MECHANICALVENTILATIONLOADSMONTHLY" ) {
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION NO LOAD HEAT REMOVAL ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION COOLING LOAD INCREASE ENERGY" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE MECHANICAL VENTILATION COOLING LOAD INCREASE DUE TO OVERHEATING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION COOLING LOAD DECREASE ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION NO LOAD HEAT ADDITION ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION HEATING LOAD INCREASE ENERGY" );
			AddRecordToOutputVariableStructure( "*",
			                                    "ZONE MECHANICAL VENTILATION HEATING LOAD INCREASE DUE TO OVERCOOLING ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION HEATING LOAD DECREASE ENERGY" );
			AddRecordToOutputVariableStructure( "*", "ZONE MECHANICAL VENTILATION AIR CHANGES PER HOUR" );
		}
	}

	void
	InputProcessor::AddRecordToOutputVariableStructure(
	std::string const & KeyValue,
	std::string const & VariableName
	) {
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   July 2010

		// PURPOSE OF THIS SUBROUTINE:
		// This routine adds a new record (if necessary) to the Output Variable
		// reporting structure.  DataOutputs, OutputVariablesForSimulation

		// METHODOLOGY EMPLOYED:
		// OutputVariablesForSimulation is a linked list structure for later
		// semi-easy perusal.

		// Using/Aliasing
		using namespace DataOutputs;

		int CurNum;
		int NextNum;
		bool FoundOne;
		std::string::size_type vnameLen; // if < length, there were units on the line/name

		std::string::size_type const rbpos = index( VariableName, '[' );
		if ( rbpos == std::string::npos ) {
			vnameLen = len_trim( VariableName );
		} else {
			vnameLen = len_trim( VariableName.substr( 0, rbpos ) );
		}

		FoundOne = false;
		std::string const VarName( VariableName.substr( 0, vnameLen ) );
		for ( CurNum = 1; CurNum <= NumConsideredOutputVariables; ++CurNum ) {
			if ( VarName == OutputVariablesForSimulation( CurNum ).VarName ) {
				FoundOne = true;
				break;
			}
		}

		if ( !FoundOne ) {
			if ( NumConsideredOutputVariables == MaxConsideredOutputVariables ) {
				ReAllocateAndPreserveOutputVariablesForSimulation();
			}
			++NumConsideredOutputVariables;
			OutputVariablesForSimulation( NumConsideredOutputVariables ).Key = KeyValue;
			OutputVariablesForSimulation( NumConsideredOutputVariables ).VarName = VarName;
			OutputVariablesForSimulation( NumConsideredOutputVariables ).Previous = 0;
			OutputVariablesForSimulation( NumConsideredOutputVariables ).Next = 0;
		} else {
			if ( KeyValue != OutputVariablesForSimulation( CurNum ).Key ) {
				NextNum = CurNum;
				if ( OutputVariablesForSimulation( NextNum ).Next != 0 ) {
					while ( OutputVariablesForSimulation( NextNum ).Next != 0 ) {
						CurNum = NextNum;
						NextNum = OutputVariablesForSimulation( NextNum ).Next;
					}
					if ( NumConsideredOutputVariables == MaxConsideredOutputVariables ) {
						ReAllocateAndPreserveOutputVariablesForSimulation();
					}
					++NumConsideredOutputVariables;
					OutputVariablesForSimulation( NumConsideredOutputVariables ).Key = KeyValue;
					OutputVariablesForSimulation( NumConsideredOutputVariables ).VarName = VarName;
					OutputVariablesForSimulation( NumConsideredOutputVariables ).Previous = NextNum;
					OutputVariablesForSimulation( NextNum ).Next = NumConsideredOutputVariables;
				} else {
					if ( NumConsideredOutputVariables == MaxConsideredOutputVariables ) {
						ReAllocateAndPreserveOutputVariablesForSimulation();
					}
					++NumConsideredOutputVariables;
					OutputVariablesForSimulation( NumConsideredOutputVariables ).Key = KeyValue;
					OutputVariablesForSimulation( NumConsideredOutputVariables ).VarName = VarName;
					OutputVariablesForSimulation( NumConsideredOutputVariables ).Previous = CurNum;
					OutputVariablesForSimulation( CurNum ).Next = NumConsideredOutputVariables;
				}
			}
		}
	}

	void
	InputProcessor::ReAllocateAndPreserveOutputVariablesForSimulation() {
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Linda Lawrie
		//       DATE WRITTEN   April 2011

		// PURPOSE OF THIS SUBROUTINE:
		// This routine does a simple reallocate for the OutputVariablesForSimulation structure, preserving
		// the data that is already in the structure.

		using namespace DataOutputs;

		int const OutputVarAllocInc( 500 );

		// up allocation by OutputVarAllocInc
		OutputVariablesForSimulation.redimension( MaxConsideredOutputVariables += OutputVarAllocInc );
	}

//	std::string
//	InputProcessor::IPTrimSigDigits( int const IntegerValue ) {
//
//		// FUNCTION INFORMATION:
//		//       AUTHOR         Linda K. Lawrie
//		//       DATE WRITTEN   March 2002
//		//       MODIFIED       na
//		//       RE-ENGINEERED  na
//
//		// PURPOSE OF THIS FUNCTION:
//		// This function accepts a number as parameter as well as the number of
//		// significant digits after the decimal point to report and returns a string
//		// that is appropriate.
//
//		std::string String; // Working string
//		gio::write( String, fmtLD ) << IntegerValue;
//		return stripped( String );
//	}
}
