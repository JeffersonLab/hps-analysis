
#include "TriggerData.h"

TriggerData::TriggerData(EVENT::LCGenericObject* trigger_data)
    : time_stamp(0), 
      single0(false), 
      single1(false),
      pair0(false),
      pair1(false),
      pulser(false) { 

    parseTriggerData(trigger_data);
}

void TriggerData::parseTriggerData(EVENT::LCGenericObject* trigger_data) { 

    int trigger_data_int = trigger_data->getIntVal(1);
    single0 = ((trigger_data_int >> 24) & 1) == 1;
    single1 = ((trigger_data_int >> 25) & 1) == 1;
    pair0 = ((trigger_data_int >> 26) & 1) == 1;
    pair1 = ((trigger_data_int >> 27) & 1) == 1;
    pulser = ((trigger_data_int >> 29) & 1) == 1;

    trigger_data_int = trigger_data->getIntVal(3);
    long w1 = trigger_data_int & 0xffffffffL;
    trigger_data_int = trigger_data->getIntVal(4);
    long w2 = trigger_data_int & 0xffffffffL;

    long timelo = w1;
    long timehi = (w2 & 0xffff) << 32;

    time_stamp = 4 * (timelo + timehi); 
}
