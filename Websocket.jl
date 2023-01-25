using HTTP.WebSockets
using JSON



server = WebSockets.listen!("127.0.0.1",2000) do ws
    for msg in ws 
        println("HI")
        send(msg,"Im back!")
    end

    
end
close(server)

