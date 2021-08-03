#include "CASourceFeature.hpp"
#include "BasicTopologyEntity.hpp"
#include "Body.hpp"
#include "RefEntityName.hpp"
#include "CastTo.hpp"
#include "CubitMessage.hpp"
#include "TDSourceFeature.hpp"

CubitAttrib* CASourceFeature_creator(RefEntity* entity, const CubitSimpleAttrib &p_csa)
{
  return new CASourceFeature(entity, p_csa);
}

CASourceFeature::CASourceFeature(RefEntity* new_attrib_owner,
                                 const CubitSimpleAttrib& csa_ptr)
                                 : CubitAttrib(new_attrib_owner)
{
    sourceFeature = GeometryFeatureEngine::FEATURE_UNDEFINED;

    if(!csa_ptr.isEmpty())
    {

      const std::vector<CubitString>& cs_list = csa_ptr.string_data_list();

      // step over the attribute type
      // now read name / option pairs
      if(cs_list.size()==2)
      {
        CubitString cs = cs_list[1];
          if (cs.length() == 0)
              PRINT_WARNING("Empty feature attribute for %s %d.\n",
              attribOwnerEntity->class_name(),
              attribOwnerEntity->id());
          else
              sourceFeature = string_to_feature_type(cs);
      }
      else
          deleteAttrib = CUBIT_TRUE;
    }
}

CASourceFeature::~CASourceFeature()
{
}

CubitStatus CASourceFeature::actuate()
{
    if (hasActuated == CUBIT_TRUE)
        return CUBIT_SUCCESS;

    // create a TDSourceFeature for the entity, if it doesn't already
    // exist
    TDSourceFeature *source_feature_data = 
        (TDSourceFeature *) attrib_owner()->get_TD(&TDSourceFeature::is_source_feature);

    if (!source_feature_data) 
    {
        // else make a new one
        TDSourceFeature* new_tool_data = new TDSourceFeature(sourceFeature);
        attrib_owner()->add_TD(new_tool_data);
    }

    delete_attrib(CUBIT_TRUE);
    hasActuated = CUBIT_TRUE;

    return CUBIT_SUCCESS;
}

CubitStatus CASourceFeature::update()
{
    if(hasUpdated)
        return CUBIT_SUCCESS;

    // set the updated flag
    hasUpdated = CUBIT_TRUE;

    // if the owner has a unique id, save it, otherwise delete this one
    TDSourceFeature *source_feature_data = 
        (TDSourceFeature *) attrib_owner()->get_TD(&TDSourceFeature::is_source_feature);

    if (!source_feature_data)
        delete_attrib(CUBIT_TRUE);
    else
    {
        if (delete_attrib() == CUBIT_TRUE)
            delete_attrib(CUBIT_FALSE);

        sourceFeature = source_feature_data->source_feature();
    }

    return CUBIT_SUCCESS;
}

CubitStatus CASourceFeature::reset()
{
    sourceFeature = GeometryFeatureEngine::FEATURE_UNDEFINED;

    hasUpdated = CUBIT_FALSE;
    return CUBIT_SUCCESS;
}

CubitSimpleAttrib CASourceFeature::cubit_simple_attrib()
{
    std::vector<CubitString> cs_list;

    // pack the string list:
    // character type of this CA
    cs_list.push_back(att_internal_name());

    // name, option pairs
    cs_list.push_back(feature_type_to_string(sourceFeature));

    return CubitSimpleAttrib(&cs_list, NULL, NULL);
}


void CASourceFeature::print()
{
    // print info on this attribute
    PRINT_INFO("CASourceFeature: owner = %s %d; feature: ",
        attribOwnerEntity->class_name(), attribOwnerEntity->id());
    PRINT_INFO("%s ", feature_type_to_string(sourceFeature).c_str());

    PRINT_INFO("\n");
}

GeometryFeatureEngine::FeatureType 
CASourceFeature::string_to_feature_type(CubitString value_in)
{
    /*
    FEATURE_UNDEFINED,
    FEATURE_HOLE,       
    FEATURE_ROUND,      
    FEATURE_CHAMFER,    
    FEATURE_SLOT ,      
    FEATURE_CUT,
    FEATURE_IMPRINT
    */
    if(value_in == "IMPRINT")
        return GeometryFeatureEngine::FEATURE_IMPRINT;
    else if(value_in == "HOLE")
        return GeometryFeatureEngine::FEATURE_HOLE;
    else if(value_in == "ROUND")
        return GeometryFeatureEngine::FEATURE_ROUND;
    else if(value_in == "CHAMFER")
        return GeometryFeatureEngine::FEATURE_CHAMFER;
    else if(value_in == "SLOT")
        return GeometryFeatureEngine::FEATURE_SLOT;
    else if(value_in == "CUT")
        return GeometryFeatureEngine::FEATURE_CUT;
    else
        return GeometryFeatureEngine::FEATURE_UNDEFINED;
}

// returning a cubit simple attribute string for the input feature type
CubitString
CASourceFeature::feature_type_to_string(GeometryFeatureEngine::FeatureType type_in)
{
    switch(type_in)
    {
    case GeometryFeatureEngine::FEATURE_IMPRINT:
        return "IMPRINT";
    case GeometryFeatureEngine::FEATURE_HOLE:
        return "HOLE";
    case GeometryFeatureEngine::FEATURE_ROUND:
        return "ROUND";
    case GeometryFeatureEngine::FEATURE_CHAMFER:
        return "CHAMFER";
    case GeometryFeatureEngine::FEATURE_SLOT:
        return "SLOT";
    case GeometryFeatureEngine::FEATURE_CUT:
        return "CUT";
    default:
        return "";
    }
}
